#include "sms.h"
#include "cache.h"

#include <algorithm>
#include <iostream>
#include <utility>
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

uint32_t SMS::prefetch(CACHE* cache, uint64_t block_number, uint32_t metadata_in)
{
  const uint64_t region            = region_of(block_number);
  const uint32_t off               = offset_of(block_number);
  const uint64_t region_base_line  = region << LOG2_REGION_BLOCKS;
  const uint64_t trigger_page      = (region_base_line << LOG2_BLOCK_SIZE) >> LOG2_PAGE_SIZE;

  auto* p = pb.find(region);
  if (!p) return 0;
  pb.touch(region);

  p->data.pending.reset(off);

  uint32_t issued = 0;
  for (uint32_t d = 1; d < REGION_BLOCKS; ++d) {
    for (int sgn = +1; sgn >= -1; sgn -= 2) {
      const int pf_off = static_cast<int>(off) + sgn * static_cast<int>(d);
      if (pf_off < 0 || pf_off >= static_cast<int>(REGION_BLOCKS)) continue;
      if (!p->data.pending.test(pf_off)) continue;

      if (cache->get_mshr_occupancy_ratio() >= MSHR_FULL_RATIO) {
        ++prefetches_blocked;
        return issued;
      }
      const auto pq_occ = cache->get_pq_occupancy_ratio();
      if (!pq_occ.empty() && pq_occ.front() >= PQ_FULL_RATIO) {
        ++prefetches_blocked;
        return issued;
      }

      const uint64_t pf_line = region_base_line + static_cast<uint64_t>(pf_off);
      const uint64_t pf_addr = pf_line << LOG2_BLOCK_SIZE;

      if ((pf_addr >> LOG2_PAGE_SIZE) != trigger_page) {
        p->data.pending.reset(pf_off);
        continue;
      }

      const bool fill_here =
        cache->get_mshr_occupancy_ratio() < MSHR_DEMOTE_RATIO;

      if (cache->prefetch_line(pf_addr, fill_here, metadata_in)) {
        p->data.pending.reset(pf_off);
        ++issued;
        ++prefetches_issued;
        if (fill_here) ++prefetches_l2_filled;
        else           ++prefetches_llc_only;
      } else {
        ++prefetches_blocked;
        return issued;
      }
    }
  }

  if (p->data.pending.none()) pb.erase(region);
  return issued;
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
  std::cout << "  Region:        " << REGION_BLOCKS << " blocks ("
            << REGION_SIZE << " B, page-aligned)" << std::endl;
  std::cout << "  Filter Table:  " << FT_ENTRIES  << " entries, "
            << FT_WAYS  << "-way" << std::endl;
  std::cout << "  Accum. Table:  " << AT_ENTRIES  << " entries, "
            << AT_WAYS  << "-way (fully associative)" << std::endl;
  std::cout << "  PHT:           " << PHT_ENTRIES << " entries, "
            << PHT_WAYS << "-way (CompactPattern: phase+period+cycle)" << std::endl;
  std::cout << "  Prefetch Buf:  " << PB_ENTRIES  << " entries, "
            << PB_WAYS  << "-way" << std::endl;
  std::cout << "  MAX_PERIOD="   << static_cast<uint32_t>(CompactPattern::MAX_PERIOD)
            << "  COVERAGE>="    << COVERAGE_MIN
            << "  SPURIOUS<="    << SPURIOUS_MAX << std::endl;
  std::cout << "  L2 fill policy: fill_here if MSHR<" << MSHR_DEMOTE_RATIO
            << ", LLC-only if " << MSHR_DEMOTE_RATIO << "<=MSHR<" << MSHR_FULL_RATIO
            << ", stop if MSHR>=" << MSHR_FULL_RATIO << std::endl;
}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip,
                                         uint8_t cache_hit, bool useful_prefetch,
                                         uint8_t type, uint32_t metadata_in)
{
  if (static_cast<size_t>(cpu) >= sms_per_cpu.size() || sms_per_cpu[cpu] == nullptr)
    return metadata_in;
  SMS* sms = sms_per_cpu[cpu];

  if (useful_prefetch) ++sms->useful_prefetches;

  const uint64_t block_number = addr >> LOG2_BLOCK_SIZE;

  sms->access(block_number, ip);
  sms->prefetch(this, block_number, metadata_in);

  return metadata_in;
}

uint32_t CACHE::prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way,
                                      uint8_t prefetch, uint64_t evicted_addr,
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

void CACHE::prefetcher_final_stats()
{
  if (static_cast<size_t>(cpu) >= sms_per_cpu.size() || sms_per_cpu[cpu] == nullptr) return;
  SMS* sms = sms_per_cpu[cpu];

  const double match_prob =
    sms->pht_lookups ? (1.0 * sms->pht_hits / sms->pht_lookups) : 0.0;
  const double avg_score =
    sms->score_count ? (sms->score_sum / sms->score_count) : 0.0;
  const double keep_rate =
    sms->comp_total ? (1.0 * sms->comp_kept / sms->comp_total) : 0.0;
  const double bit_coverage =
    sms->set_bits_total
      ? (1.0 * sms->set_bits_reconstructed / sms->set_bits_total) : 0.0;
  const double bit_coverage_of_kept =
    sms->set_bits_kept
      ? (1.0 * sms->set_bits_reconstructed / sms->set_bits_kept) : 0.0;
  const double bit_spur_rate =
    sms->set_bits_kept
      ? (1.0 * sms->set_bits_spurious / sms->set_bits_kept) : 0.0;
  const double compression =
    sms->storage_bits_actual
      ? (1.0 * sms->storage_bits_raw_basis / sms->storage_bits_actual) : 0.0;

  const double pf_useful_rate =
    sms->prefetches_issued
      ? (1.0 * sms->useful_prefetches / sms->prefetches_issued) : 0.0;

  std::cout << NAME << " SMS-Fourier stats:"
            << " pht_lookups="  << sms->pht_lookups
            << " pht_hits="     << sms->pht_hits
            << " match_prob="   << match_prob
            << " generations="  << sms->generations
            << " pf_issued="    << sms->prefetches_issued
            << " pf_l2_fill="   << sms->prefetches_l2_filled
            << " pf_llc_only="  << sms->prefetches_llc_only
            << " pf_blocked="   << sms->prefetches_blocked
            << " pf_useful="    << sms->useful_prefetches
            << " pf_useful_r="  << pf_useful_rate
            << std::endl;

  std::cout << NAME << " SMS-Fourier compression:"
            << " total="            << sms->comp_total
            << " kept="             << sms->comp_kept
            << " keep_rate="        << keep_rate
            << " drop_empty="       << sms->comp_dropped_empty
            << " drop_no_period="   << sms->comp_dropped_no_period
            << " drop_period_big="  << sms->comp_dropped_period_large
            << " drop_coverage="    << sms->comp_dropped_coverage
            << " avg_score="        << avg_score
            << std::endl;

  std::cout << NAME << " SMS-Fourier period_hist:";
  for (std::size_t p = 1; p <= CompactPattern::MAX_PERIOD; ++p) {
    std::cout << " p" << p << "=" << sms->comp_period_hist[p];
  }
  std::cout << std::endl;

  std::cout << NAME << " SMS-Fourier coverage:"
            << " set_bits_total="    << sms->set_bits_total
            << " set_bits_kept="     << sms->set_bits_kept
            << " reconstructed="     << sms->set_bits_reconstructed
            << " spurious="          << sms->set_bits_spurious
            << " bit_cov_overall="   << bit_coverage
            << " bit_cov_of_kept="   << bit_coverage_of_kept
            << " bit_spur_rate="     << bit_spur_rate
            << std::endl;

  std::cout << NAME << " SMS-Fourier storage:"
            << " raw_basis_bits=" << sms->storage_bits_raw_basis
            << " actual_bits="    << sms->storage_bits_actual
            << " compression_x="  << compression
            << std::endl;

  static const char* const POP_LABELS[SMS::POP_BUCKETS] = {
    "0", "1", "2", "3-4", "5-8", "9-16", "17-32", "33-64"
  };
  std::cout << NAME << " SMS-Fourier popcount_hist:";
  for (std::size_t b = 0; b < SMS::POP_BUCKETS; ++b) {
    std::cout << " " << POP_LABELS[b] << "=" << sms->popcount_hist[b];
  }
  std::cout << std::endl;

  constexpr std::size_t TOP_N = 15;
  std::vector<std::pair<uint16_t, uint64_t>> ranked(
      sms->stencil_hist.begin(), sms->stencil_hist.end());
  std::sort(ranked.begin(), ranked.end(),
            [](const std::pair<uint16_t, uint64_t>& a,
               const std::pair<uint16_t, uint64_t>& b) {
              if (a.second != b.second) return a.second > b.second;
              return a.first < b.first;
            });
  const std::size_t n_show = std::min(TOP_N, ranked.size());
  std::cout << NAME << " SMS-Fourier top_stencils"
            << " (unique=" << ranked.size() << ", showing top " << n_show << "):";
  for (std::size_t i = 0; i < n_show; ++i) {
    const uint8_t period = static_cast<uint8_t>(ranked[i].first >> 8);
    const uint8_t cycle  = static_cast<uint8_t>(ranked[i].first & 0xFF);
    std::cout << " p" << static_cast<uint32_t>(period) << "/";
    for (uint8_t k = 0; k < period; ++k) {
      std::cout << static_cast<uint32_t>((cycle >> k) & 1u);
    }
    std::cout << "=" << ranked[i].second;
  }
  std::cout << std::endl;
}
