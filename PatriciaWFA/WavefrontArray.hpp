#pragma once

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <numeric>

class WavefrontArray {
 public:
    // =========================================================
    // Constants for vk encoding
    // =========================================================
    static constexpr uint64_t DIAGONAL_MASK = 0xFFFFFFFF;
    static constexpr int32_t DIAGONAL_OFFSET = 0x40000000;

    struct DiagState {
        uint64_t vk;
        int32_t offset;
    };

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
    
    // Pre-allocates capacity for the internal vectors (alternative to vector::reserve).
    // Used to prevent unintended reallocations (realloc) during loops as the wavefront expands.
    void reserve_capacity(size_t expected_elements) {
        _states.reserve(expected_elements);
    }

    // Resets the logical size to 0 while preserving the allocated physical memory capacity (alternative to vector::clear).
    // Expected to be called when the instance is reused and retrieved from the WavefrontPool.
    void clear_logical_size() {
        _active_size = 0;
    }

    void set_size(size_t new_size) {
        _active_size = new_size;
        if (_states.size() < new_size) {
            // This is not preferred usage
            _states.resize(new_size);
        }
    }

    // =========================================================
    // 2. Element Addition & Scalar Access
    // =========================================================

    // 読み取り（値返しによる安全性の担保）
    uint64_t get_vk(size_t idx) const noexcept { return _states[idx].vk; }
    int32_t get_offset(size_t idx) const noexcept { return _states[idx].offset; }

    // 書き込み（エンコードロジックの強制）
    void update_state(size_t idx, uint32_t node_id, int32_t k, int32_t offset) noexcept {
        _states[idx] = {calc_vk(node_id, k), offset};
    }

    void update_state(size_t idx, uint64_t vk, int32_t offset) noexcept {
        _states[idx] = {vk, offset};
    }

    // 新規状態の追加 (node_id と k から vk を計算)
    void push_back_state(uint32_t node_id, int32_t k, int32_t offset) {
        push_back_state(calc_vk(node_id, k), offset);
    }

    // 既存の vk を直接追加するオーバーロード (マージ・フラッシュ処理での再エンコード防止用)
    void push_back_state(uint64_t vk, int32_t offset) {
        // if (_active_size > 0 && _states[_active_size - 1].vk == vk) {
        //     if (offset > _states[_active_size - 1].offset) {
        //         _states[_active_size - 1].offset = offset;
        //     }
        //     std::cout << "[WARN] Duplicate state detected in WavefrontArray::push_back_state. Merging offsets.\n";
        //     return;
        // }
        if (_active_size == _states.size()) { [[unlikely]]
            _states.push_back({vk, offset});
        } else { [[likely]]
            _states[_active_size] = {vk, offset};
        }
        _active_size++;
    }

    // Retrieves the current number of active elements (alternative to vector::size).
    size_t active_size() const {
        return _active_size;
    }

    bool empty() const {
        return _active_size == 0;
    }

    void swap(WavefrontArray& other) noexcept {
        _states.swap(other._states);
        std::swap(_active_size, other._active_size);
    }

    // =========================================================
    // 5. Utilities (Class Methods)
    // =========================================================
    
    // Restores node_id from vk (extracts the upper 32 bits).
    static uint32_t calc_node_id_from_vk(uint64_t vk) {
        return static_cast<uint32_t>(vk >> 32);
    }
    
    // Restores the diagonal k from vk (applies DIAGONAL_MASK and subtracts DIAGONAL_OFFSET).
    static int32_t calc_k_from_vk(uint64_t vk) {
        return static_cast<int32_t>(vk & DIAGONAL_MASK) - DIAGONAL_OFFSET;
    }

    // Synthesizes vk from node_id and diagonal k (adds DIAGONAL_OFFSET and performs bit shifting).
    static uint64_t calc_vk(uint32_t node_id, int32_t k) {
        return (static_cast<uint64_t>(node_id) << 32) | (static_cast<uint64_t>(k + DIAGONAL_OFFSET) & DIAGONAL_MASK);
    }

 private:
    std::vector<DiagState> _states;
    size_t _active_size = 0;
};
