#pragma once

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace dt_patricia::internal {

class WavefrontArray {
 public:
    // =========================================================
    // Constants for vk encoding
    // =========================================================
    static constexpr uint64_t DIAGONAL_MASK = 0xFFFFFFFF;
    static constexpr int32_t DIAGONAL_OFFSET = 0x40000000;

    // =========================================================
    // 0. Rule of Five
    // =========================================================
    WavefrontArray() = default;
    ~WavefrontArray() = default;
    WavefrontArray(const WavefrontArray&) = delete;
    WavefrontArray& operator=(const WavefrontArray&) = delete;
    WavefrontArray(WavefrontArray&&) noexcept = default;
    WavefrontArray& operator=(WavefrontArray&&) noexcept = default;

    // =========================================================
    // 1. Memory Management & Lifecycle
    // =========================================================

    // Pre-allocates capacity for the internal vectors.
    // Uses resize (not reserve) so physical size tracks what we can safely write to.
    void reserve_capacity(size_t n) {
        if (n > _vks.size()) {
            _vks.resize(n);
            _offsets.resize(n);
        }
    }

    // Resets the logical size to 0 while preserving the allocated physical memory capacity.
    void clear_logical_size() {
        _active_size = 0;
    }

    void set_size(size_t new_size) {
        _active_size = new_size;
        if (_vks.size() < new_size) {
            // This is not preferred usage
            _vks.resize(new_size);
            _offsets.resize(new_size);
        }
    }

    // =========================================================
    // 2. Element Addition & Scalar Access
    // =========================================================

    uint64_t get_vk(size_t idx) const noexcept { return _vks[idx]; }
    int32_t get_offset(size_t idx) const noexcept { return _offsets[idx]; }

    void update_state(size_t idx, uint32_t node_id, int32_t k, int32_t offset) noexcept {
        _vks[idx] = calc_vk(node_id, k);
        _offsets[idx] = offset;
    }

    void update_state(size_t idx, uint64_t vk, int32_t offset) noexcept {
        _vks[idx] = vk;
        _offsets[idx] = offset;
    }

    void push_back_state(uint32_t node_id, int32_t k, int32_t offset) {
        push_back_state(calc_vk(node_id, k), offset);
    }

    void push_back_state(uint64_t vk, int32_t offset) {
        if (_active_size == _vks.size()) { [[unlikely]]
            _vks.push_back(vk);
            _offsets.push_back(offset);
        } else { [[likely]]
            _vks[_active_size] = vk;
            _offsets[_active_size] = offset;
        }
        _active_size++;
    }

    // Fast path: safe to call ONLY when _active_size < _vks.size()
    // (enforced by prior reserve_capacity call).
    void push_back_unchecked(uint64_t vk, int32_t offset) noexcept {
        assert(_active_size < _vks.size());
        _vks[_active_size] = vk;
        _offsets[_active_size] = offset;
        _active_size++;
    }

    size_t active_size() const {
        return _active_size;
    }

    bool empty() const {
        return _active_size == 0;
    }

    void swap(WavefrontArray& other) noexcept {
        _vks.swap(other._vks);
        _offsets.swap(other._offsets);
        std::swap(_active_size, other._active_size);
    }

    // =========================================================
    // 5. Utilities (Class Methods)
    // =========================================================

    static uint32_t calc_node_id_from_vk(uint64_t vk) {
        return static_cast<uint32_t>(vk >> 32);
    }

    static int32_t calc_k_from_vk(uint64_t vk) {
        return static_cast<int32_t>(vk & DIAGONAL_MASK) - DIAGONAL_OFFSET;
    }

    static uint64_t calc_vk(uint32_t node_id, int32_t k) {
        return (static_cast<uint64_t>(node_id) << 32) | (static_cast<uint64_t>(k + DIAGONAL_OFFSET) & DIAGONAL_MASK);
    }

 private:
    std::vector<uint64_t> _vks;
    std::vector<int32_t> _offsets;
    size_t _active_size = 0;
};

}  // namespace dt_patricia::internal
