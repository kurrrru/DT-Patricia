#pragma once

#include <cstdint>
#include <stdexcept>

namespace dt_patricia {

struct UnitCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit = true;
    static constexpr uint32_t mismatch = 1;
    static constexpr uint32_t gap = 1;
};

struct LinearGapCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit = false;
    uint32_t mismatch;
    uint32_t gap;

    LinearGapCost(uint32_t m, uint32_t g) : mismatch(m), gap(g) {
        if (m < 1 || g < 1) {
            throw std::invalid_argument("LinearGapCost: mismatch and gap must be >= 1");
        }
    }
};

// cost(L) = gap_open + gap_extend * L
struct AffineGapCost {
    static constexpr bool is_linear = false;
    static constexpr bool is_unit = false; 
    uint32_t mismatch;
    uint32_t gap_open;
    uint32_t gap_extend;

    AffineGapCost(uint32_t m, uint32_t g_open, uint32_t g_extend)
        : mismatch(m), gap_open(g_open), gap_extend(g_extend) {
        if (m < 1 || g_extend < 1) {
            throw std::invalid_argument("AffineGapCost: mismatch and gap_extend must be >= 1");
        }
    }
};

}  // namespace dt_patricia
