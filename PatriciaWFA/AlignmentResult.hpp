#pragma once

#include <cstdint>

struct AlignmentResult {
    uint32_t string_id;  // ID of the matched string in the original data
    uint32_t score;      // Edit distance score
};
