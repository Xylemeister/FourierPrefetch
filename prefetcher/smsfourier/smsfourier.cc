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

CompactPattern SMS::compress(const SpatialPattern& /*pat*/, uint8_t /*trigger_offset*/)
{
  return CompactPattern{};
}

SpatialPattern SMS::find_in_pht(uint64_t /*pc*/, uint32_t /*offset*/)
{
  return SpatialPattern{};
}

void SMS::insert_in_pht(uint64_t /*pc*/, uint32_t /*offset*/, const SpatialPattern& /*pat*/)
{
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
                                         uint8_t /*cache_hit*/, bool /*useful_prefetch*/,
                                         uint8_t /*type*/, uint32_t metadata_in)
{
  if (static_cast<size_t>(cpu) >= sms_per_cpu.size() || sms_per_cpu[cpu] == nullptr)
    return metadata_in;
  SMS* sms = sms_per_cpu[cpu];

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
