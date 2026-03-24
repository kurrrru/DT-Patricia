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

private:
    // SoA (Structure of Arrays) 
    std::vector<uint64_t> _vks;
    std::vector<int32_t> _offsets;
    
    // Logical size of the wavefront array (number of active states). 
    // This is not necessarily equal to the physical size (capacity) of the vectors.
    size_t _active_size = 0;

public:
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
        _vks.reserve(expected_elements);
        _offsets.reserve(expected_elements);
    }

    // Resets the logical size to 0 while preserving the allocated physical memory capacity (alternative to vector::clear).
    // Expected to be called when the instance is reused and retrieved from the WavefrontPool.
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
    
    // Appends a wavefront element to the end of the arrays and increments the logical size (alternative to vector::push_back).
    // Internally writes values to both _vks and _offsets.
    void push_back_state(uint32_t node_id, int32_t k, int32_t offset) {
        const uint64_t vk = calc_vk(node_id, k);

        if (_active_size > 0 && _vks[_active_size - 1] == vk) {//[NOTE]ここのif文まるごと削除しても良さそう。後で確認
            std::cout << "Dump" << std::endl;
            for (size_t i = 0; i < _active_size; ++i) {
                uint32_t existing_node_id = calc_node_id_from_vk(_vks[i]);
                int32_t existing_k = calc_k_from_vk(_vks[i]);
                int32_t existing_offset = _offsets[i];
                std::cout << "  State " << i << ": (node_id=" << existing_node_id 
                          << ", k=" << existing_k << ", offset=" << existing_offset << ")\n";
            }
            std::cout << "Attempting to add duplicate state: (node_id=" << node_id << ", k=" << k << ", offset=" << offset << ")\n";
            if (offset > _offsets[_active_size - 1]) {
                _offsets[_active_size - 1] = offset;
            }
            // This warning is for debugging purposes to detect unintended duplicate insertions.
            // This line will be removed in the future
            
            std::cout << "[WARN] Duplicate state detected in WavefrontArray::push_back_state. Merging offsets.\n";
            return;
        }

        // 既存の確保済みキャパシティの再利用、または新規追加
        if (_active_size < _vks.size()) {
            _vks[_active_size] = vk;
            _offsets[_active_size] = offset;
        } else {
            _vks.push_back(vk);
            _offsets.push_back(offset);
        }
        
        ++_active_size;
    }

    void push_back_state(uint64_t vk, int32_t offset) {
        if (_active_size > 0 && _vks[_active_size - 1] == vk) {// [NOTE]ここも確認
            if (offset > _offsets[_active_size - 1]) {
                _offsets[_active_size - 1] = offset;
            }
            std::cout << "[WARN] Duplicate state detected in WavefrontArray::push_back_state2. Merging offsets.\n";
            return;
        }

        if (_active_size < _vks.size()) {
            _vks[_active_size] = vk;
            _offsets[_active_size] = offset;
        } else {
            _vks.push_back(vk);
            _offsets.push_back(offset);
        }
        
        ++_active_size;
    }

    // Retrieves the current number of active elements (alternative to vector::size).
    size_t active_size() const {
        return _active_size;
    }

    bool empty() const {
        return _active_size == 0;
    }

    // Retrieves the vk at the specified index (for scalar access).
    uint64_t get_vk(size_t index) const {
        return _vks[index];
    }

    // Retrieves the offset at the specified index (for scalar access).
    int32_t get_offset(size_t index) const {
        return _offsets[index];
    }

        // Swaps the contents of this WavefrontArray with another instance. This is used to efficiently replace the current wavefront array with a new one without copying elements.
    void swap(WavefrontArray& other) {
        std::swap(_vks, other._vks);
        std::swap(_offsets, other._offsets);
        std::swap(_active_size, other._active_size);
    }

    void update_state(size_t index, uint32_t node_id, int32_t k, int32_t offset) {
        const uint64_t vk = calc_vk(node_id, k);
        _vks[index] = vk;
        _offsets[index] = offset;
    }


    // =========================================================
    // 3. Pointer Access for Vector/SIMD Processing
    // =========================================================
    
    // Retrieves the start pointer of the vk array. For SIMD loads and contiguous memory access (alternative to vector::data).
    const uint64_t* data_vks() const {
        return _vks.data();
    }

    // Retrieves the start pointer of the offset array.
    const int32_t* data_offsets() const {
        return _offsets.data();
    }

    // Retrieves the mutable pointer of the offset array. For SIMD stores and batch addition processing in the Compute phase.
    int32_t* data_offsets_mutable() {
        return _offsets.data();
    }

    // =========================================================
    // 4. Algorithmic Operations
    // =========================================================
    
    // Sorts the elements within the range of the logical size (_active_size) using _vks as the primary key and eliminates duplicates.
    // * Mandatory requirement: When swapping or moving elements in _vks, the elements in _offsets at the corresponding indices must be moved synchronously.
    void sort_and_deduplicate();
    
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
};

inline void WavefrontArray::sort_and_deduplicate() {
    if (_active_size <= 1) return;

    // 1. Create an index array [0, 1, 2, ..., _active_size - 1]
    std::vector<uint32_t> indices(_active_size);
    std::iota(indices.begin(), indices.end(), 0);

    // 2. Sort indices based on the corresponding vk values in _vks
    std::sort(indices.begin(), indices.end(),
        [this](uint32_t a, uint32_t b) {
            return _vks[a] < _vks[b];
        });

    // 3. In-place permutation
    for (uint32_t i = 0; i < _active_size; ++i) {
        if (indices[i] != i) {
            uint64_t current_vk = _vks[i];
            int32_t current_offset = _offsets[i];
            uint32_t j = i;

            while (indices[j] != i) {
                uint32_t next = indices[j];
                _vks[j] = _vks[next];
                _offsets[j] = _offsets[next];
                indices[j] = j; // 処理済みマーク
                j = next;
            }
            _vks[j] = current_vk;
            _offsets[j] = current_offset;
            indices[j] = j;
        }
    }

    // 4. Deduplication
    size_t unique_size = 1;
    for (size_t i = 1; i < _active_size; ++i) {
        if (_vks[i] != _vks[unique_size - 1]) {
            _vks[unique_size] = _vks[i];
            _offsets[unique_size] = _offsets[i];
            unique_size++;
        } else {
            if (_offsets[i] > _offsets[unique_size - 1]) {
                _offsets[unique_size - 1] = _offsets[i];
            }
        }
    }
    _active_size = unique_size;
}
