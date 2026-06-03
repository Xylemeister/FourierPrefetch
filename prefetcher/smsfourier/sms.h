#ifndef SMS_FOURIER_H_
#define SMS_FOURIER_H_

#include "cache.h"
#include "transform.h"

#include <bitset>
#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace smsf_space
{
  constexpr uint32_t REGION_SIZE         = 4 * 1024;
  constexpr uint32_t REGION_BLOCKS       = REGION_SIZE / 64;
  constexpr uint32_t LOG2_REGION_BLOCKS  = 6;

  constexpr uint32_t PC_WIDTH            = 16;
  constexpr uint32_t ADDR_WIDTH          = LOG2_REGION_BLOCKS;

  constexpr uint32_t PHT_ENTRIES         = 1024;
  constexpr uint32_t PHT_WAYS            = 4;
  constexpr uint32_t FT_ENTRIES          = 64;
  constexpr uint32_t FT_WAYS             = 8;
  constexpr uint32_t AT_ENTRIES          = 64;
  constexpr uint32_t AT_WAYS             = 64;
  constexpr uint32_t PB_ENTRIES          = 32;
  constexpr uint32_t PB_WAYS             = 8;

  constexpr double   PQ_FULL_RATIO       = 0.75;
  constexpr double   MSHR_FULL_RATIO     = 0.75;
  constexpr double   MSHR_DEMOTE_RATIO   = 0.5;

  constexpr double COVERAGE_MIN = 0.80;
  constexpr double SPURIOUS_MAX = 0.25;

  constexpr bool PHASE_FROM_TRIGGER_OFFSET = false;

  using SpatialPattern = std::bitset<REGION_BLOCKS>;

  struct CompactPattern {
    static constexpr uint8_t MAX_PERIOD = 8;

    bool     valid  = false;
    uint8_t  phase  = 0;
    uint8_t  period = 0;
    uint8_t  cycle  = 0;

    SpatialPattern expand() const;

    std::size_t storage_bits() const {
      if (!valid) return 0;
      return 1 + 6 + 3 + std::size_t{period};
    }
  };

  template <typename Data>
  class SetAssocCache {
  public:
    struct Entry {
      bool     valid = false;
      uint64_t key   = 0;
      Data     data{};
      uint64_t lru   = 0;
    };

    SetAssocCache(uint32_t n_entries, uint32_t n_ways)
      : num_ways(n_ways),
        num_sets(n_entries / n_ways),
        table(num_sets, std::vector<Entry>(n_ways))
    {
      assert(n_entries % n_ways == 0);
      assert(num_sets > 0);
    }

    Entry* find(uint64_t key)
    {
      for (auto& e : table[set_of(key)]) {
        if (e.valid && e.key == key) return &e;
      }
      return nullptr;
    }

    Entry insert(uint64_t key, const Data& data)
    {
      auto& set = table[set_of(key)];
      Entry* victim     = nullptr;
      uint64_t min_lru  = UINT64_MAX;

      for (auto& e : set) {
        if (e.valid && e.key == key) {
          e.data = data;
          e.lru  = ++lru_counter;
          return Entry{};
        }
        if (!e.valid) {
          e.valid = true;
          e.key   = key;
          e.data  = data;
          e.lru   = ++lru_counter;
          return Entry{};
        }
        if (e.lru < min_lru) { min_lru = e.lru; victim = &e; }
      }
      Entry evicted = *victim;
      victim->key  = key;
      victim->data = data;
      victim->lru  = ++lru_counter;
      return evicted;
    }

    void touch(uint64_t key)
    {
      if (Entry* e = find(key)) e->lru = ++lru_counter;
    }

    Entry erase(uint64_t key)
    {
      for (auto& e : table[set_of(key)]) {
        if (e.valid && e.key == key) {
          Entry erased = e;
          e.valid = false;
          return erased;
        }
      }
      return Entry{};
    }

  private:
    uint32_t set_of(uint64_t key) const
    {
      uint64_t h = key;
      h ^= h >> 21;
      h ^= h >> 13;
      h ^= h >> 7;
      return static_cast<uint32_t>(h % num_sets);
    }

    uint32_t num_ways;
    uint32_t num_sets;
    std::vector<std::vector<Entry>> table;
    uint64_t lru_counter = 0;
  };

  struct FTData {
    uint64_t pc;
    uint32_t offset;
  };

  struct ATData {
    uint64_t       pc;
    uint32_t       offset;
    SpatialPattern pattern;
  };

  struct PHTData {
    CompactPattern pattern;
  };

  struct PBData {
    SpatialPattern pending;
  };

  class SMS {
  public:
    SMS();

    void access(uint64_t block_number, uint64_t pc);
    uint32_t prefetch(CACHE* cache, uint64_t block_number, uint32_t metadata_in);
    void eviction(uint64_t block_number);

    uint64_t pht_lookups        = 0;
    uint64_t pht_hits           = 0;
    uint64_t generations        = 0;
    uint64_t prefetches_issued    = 0;
    uint64_t prefetches_blocked   = 0;
    uint64_t prefetches_l2_filled = 0;
    uint64_t prefetches_llc_only  = 0;
    uint64_t useful_prefetches    = 0;

    uint64_t comp_total                = 0;
    uint64_t comp_kept                 = 0;
    uint64_t comp_dropped_empty        = 0;
    uint64_t comp_dropped_no_period    = 0;
    uint64_t comp_dropped_period_large = 0;
    uint64_t comp_dropped_coverage     = 0;

    uint64_t comp_period_hist[CompactPattern::MAX_PERIOD + 1] = {0};

    std::unordered_map<uint16_t, uint64_t> stencil_hist;

    static constexpr std::size_t POP_BUCKETS = 8;
    uint64_t popcount_hist[POP_BUCKETS] = {0};

    uint64_t set_bits_total           = 0;
    uint64_t set_bits_kept            = 0;
    uint64_t set_bits_reconstructed   = 0;
    uint64_t set_bits_spurious        = 0;

    uint64_t storage_bits_raw_basis   = 0;
    uint64_t storage_bits_actual      = 0;

    double   score_sum                = 0.0;
    uint64_t score_count              = 0;

  private:
    SetAssocCache<FTData>  ft;
    SetAssocCache<ATData>  at;
    SetAssocCache<PHTData> pht;
    SetAssocCache<PBData>  pb;

    SpatialPattern find_in_pht(uint64_t pc, uint32_t offset);
    void           insert_in_pht(uint64_t pc, uint32_t offset, const SpatialPattern& pat);

    CompactPattern compress(const SpatialPattern& pat, uint8_t trigger_offset);

    static uint64_t region_of(uint64_t block_number) {
      return block_number >> LOG2_REGION_BLOCKS;
    }
    static uint32_t offset_of(uint64_t block_number) {
      return static_cast<uint32_t>(block_number & (REGION_BLOCKS - 1));
    }
    static uint64_t make_pht_key(uint64_t pc, uint32_t offset) {
      pc     &= (1ULL << PC_WIDTH)   - 1;
      offset &= (1U   << ADDR_WIDTH) - 1;
      return (pc << ADDR_WIDTH) | offset;
    }
  };
} // namespace smsf_space

#endif // SMS_FOURIER_H_
