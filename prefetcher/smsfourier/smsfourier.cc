#include "sms.h"
#include "cache.h"

#include <iostream>
#include <vector>

using namespace smsf_space;

namespace {
  std::vector<SMS*> sms_per_cpu;
}

SpatialPattern CompactPattern::expand() const
{
  SpatialPattern out;
  if (!valid || period == 0) return out;
  for (uint32_t i = phase; i < REGION_BLOCKS; ++i) {
    const uint32_t k = (i - phase) % period;
    if ((cycle >> k) & 1u) out.set(i);
  }
  return out;
}

SMS::SMS()
  : ft (FT_ENTRIES,  FT_WAYS),
    at (AT_ENTRIES,  AT_WAYS),
    pht(PHT_ENTRIES, PHT_WAYS),
    pb (PB_ENTRIES,  PB_WAYS)
{}

CompactPattern SMS::compress(const SpatialPattern& pat, uint8_t trigger_offset)
{
  CompactPattern out;
  ++comp_total;
  storage_bits_raw_basis += REGION_BLOCKS;

  const std::size_t pop = pat.count();
  set_bits_total += pop;

  {
    std::size_t b;
    if      (pop == 0)  b = 0;
    else if (pop == 1)  b = 1;
    else if (pop == 2)  b = 2;
    else if (pop <= 4)  b = 3;
    else if (pop <= 8)  b = 4;
    else if (pop <= 16) b = 5;
    else if (pop <= 32) b = 6;
    else                b = 7;
    ++popcount_hist[b];
  }

  if (pop == 0) {
    ++comp_dropped_empty;
    return out;
  }

  periodicity::PeriodicityBuf<REGION_BLOCKS> buf;
  for (uint32_t i = 0; i < REGION_BLOCKS; ++i) {
    buf.Insert(pat.test(i) ? 1 : 0);
  }

  auto pred = buf.Predict(
      periodicity::PERIODICITY_THRESHOLD,
      periodicity::PERIODICITY_L2_THRESHOLD,
      1,
      0,
      0);

  if (pred.period == 0) {
    ++comp_dropped_no_period;
    return out;
  }
  if (pred.period > CompactPattern::MAX_PERIOD) {
    ++comp_dropped_period_large;
    return out;
  }

  uint8_t phase = 0;
  if constexpr (PHASE_FROM_TRIGGER_OFFSET) {
    phase = trigger_offset;
  } else {
    for (uint32_t i = 0; i < REGION_BLOCKS; ++i) {
      if (pat.test(i)) { phase = static_cast<uint8_t>(i); break; }
    }
  }

  uint8_t cycle = 0;
  for (uint32_t k = 0; k < pred.period; ++k) {
    const uint32_t pos = static_cast<uint32_t>(phase) + k;
    if (pos < REGION_BLOCKS && pat.test(pos)) cycle |= static_cast<uint8_t>(1u << k);
  }

  CompactPattern cand;
  cand.valid  = true;
  cand.phase  = phase;
  cand.period = static_cast<uint8_t>(pred.period);
  cand.cycle  = cycle;

  const SpatialPattern rec  = cand.expand();
  const std::size_t   hit  = (rec & pat).count();
  const std::size_t   spur = (rec & ~pat).count();

  ++comp_kept;
  ++comp_period_hist[cand.period];
  const uint16_t stencil_key =
      static_cast<uint16_t>((static_cast<uint16_t>(cand.period) << 8) | cand.cycle);
  ++stencil_hist[stencil_key];
  set_bits_kept           += pop;
  set_bits_reconstructed  += hit;
  set_bits_spurious       += spur;
  storage_bits_actual     += cand.storage_bits();
  score_sum               += pred.best_score;
  ++score_count;

  return cand;
}

SpatialPattern SMS::find_in_pht(uint64_t pc, uint32_t offset)
{
  ++pht_lookups;
  uint64_t key = make_pht_key(pc, offset);
  if (auto* e = pht.find(key)) {
    pht.touch(key);
    ++pht_hits;
    return e->data.pattern.expand();
  }
  return SpatialPattern{};
}

void SMS::insert_in_pht(uint64_t pc, uint32_t offset, const SpatialPattern& pat)
{
  CompactPattern compact = compress(pat, static_cast<uint8_t>(offset));
  if (!compact.valid) return;
  pht.insert(make_pht_key(pc, offset), PHTData{compact});
}

void SMS::access(uint64_t block_number, uint64_t pc)
{
  const uint64_t region = region_of(block_number);
  const uint32_t off    = offset_of(block_number);

  if (auto* a = at.find(region)) {
    a->data.pattern.set(off);
    at.touch(region);
    return;
  }

  auto* f = ft.find(region);
  if (!f) {
    ft.insert(region, FTData{pc, off});

    SpatialPattern pred = find_in_pht(pc, off);
    if (pred.any()) {
      pred.reset(off);
      pb.insert(region, PBData{pred});
    }
    return;
  }

  if (f->data.offset != off) {
    ATData ae;
    ae.pc     = f->data.pc;
    ae.offset = f->data.offset;
    ae.pattern.set(f->data.offset);
    ae.pattern.set(off);

    auto victim = at.insert(region, ae);
    ft.erase(region);

    if (victim.valid) {
      insert_in_pht(victim.data.pc, victim.data.offset, victim.data.pattern);
      ++generations;
    }
  }
}

uint32_t SMS::prefetch(CACHE* /*cache*/, uint64_t /*block_number*/, uint32_t /*metadata_in*/)
{
  return 0;
}

void SMS::eviction(uint64_t block_number)
{
  const uint64_t region = region_of(block_number);
  ft.erase(region);

  auto victim = at.erase(region);
  if (victim.valid) {
    insert_in_pht(victim.data.pc, victim.data.offset, victim.data.pattern);
    ++generations;
  }
}

void CACHE::prefetcher_initialize()
{
  if (sms_per_cpu.size() <= static_cast<size_t>(cpu))
    sms_per_cpu.resize(static_cast<size_t>(cpu) + 1, nullptr);
  sms_per_cpu[cpu] = new SMS();

  std::cout << NAME << " SMS-Fourier Prefetcher (L2)" << std::endl;
}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip,
                                         uint8_t /*cache_hit*/, bool useful_prefetch,
                                         uint8_t /*type*/, uint32_t metadata_in)
{
  if (static_cast<size_t>(cpu) >= sms_per_cpu.size() || sms_per_cpu[cpu] == nullptr)
    return metadata_in;
  SMS* sms = sms_per_cpu[cpu];

  if (useful_prefetch) ++sms->useful_prefetches;

  const uint64_t block_number = addr >> LOG2_BLOCK_SIZE;
  sms->access(block_number, ip);

  return metadata_in;
}

uint32_t CACHE::prefetcher_cache_fill(uint64_t /*addr*/, uint32_t /*set*/, uint32_t /*way*/,
                                      uint8_t /*prefetch*/, uint64_t evicted_addr,
                                      uint32_t metadata_in)
{
  if (static_cast<size_t>(cpu) >= sms_per_cpu.size() || sms_per_cpu[cpu] == nullptr)
    return metadata_in;
  SMS* sms = sms_per_cpu[cpu];

  if (evicted_addr != 0) {
    sms->eviction(evicted_addr >> LOG2_BLOCK_SIZE);
  }
  return metadata_in;
}

void CACHE::prefetcher_cycle_operate() {}

void CACHE::prefetcher_final_stats() {}
