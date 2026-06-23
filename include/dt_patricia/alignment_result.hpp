#pragma once

#include <cstdint>

namespace dt_patricia {

struct AlignmentResult {
    uint32_t string_id;  // ID of the matched string in the original data
    uint32_t score;      // Edit distance score
};

}  // namespace dt_patricia
