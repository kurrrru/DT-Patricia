#pragma once

#include <cstdint>

namespace dt_patricia::internal {

// 常に a < mod, b < mod であることを前提とする高速な剰余加算
[[nodiscard]] inline uint32_t add_mod(uint32_t a, uint32_t b, uint32_t mod) noexcept {
    uint32_t res = a + b;
    return (res >= mod) ? res - mod : res;
}

// 常に a < mod, b < mod であることを前提とする高速な剰余減算
[[nodiscard]] inline uint32_t sub_mod(uint32_t a, uint32_t b, uint32_t mod) noexcept {
    return (a >= b) ? a - b : a + mod - b;
}

[[nodiscard]] inline uint32_t increment_mod(uint32_t a, uint32_t mod) noexcept {
    return (a + 1 >= mod) ? 0 : a + 1;
}

[[nodiscard]] inline uint32_t decrement_mod(uint32_t a, uint32_t mod) noexcept {
    return (a == 0) ? mod - 1 : a - 1;
}

}  // namespace dt_patricia::internal
