#ifndef SMS_FOURIER_H_
#define SMS_FOURIER_H_

#include "cache.h"
#include "transform.h"

#include <bitset>
#include <cstddef>
#include <cstdint>

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
} // namespace smsf_space

#endif // SMS_FOURIER_H_
