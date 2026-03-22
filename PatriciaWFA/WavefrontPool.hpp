#pragma once

#include "WavefrontArray.hpp"
#include <vector>
#include <memory>
#include <cstddef>

class WavefrontPool {
private:
    std::vector<std::unique_ptr<WavefrontArray>> _pool;
    
    // Index of the next available element in the pool
    size_t _cursor = 0;

public:
    // =========================================================
    // 0. Rule of Five
    // =========================================================
    
    explicit WavefrontPool(size_t initial_pool_size = 128) {
        _pool.reserve(initial_pool_size);
        for (size_t i = 0; i < initial_pool_size; ++i) {
            _pool.push_back(std::make_unique<WavefrontArray>());
        }
    }
    
    ~WavefrontPool() = default;
    
    // Prevent unintended copying of the pool itself
    WavefrontPool(const WavefrontPool&) = delete;
    WavefrontPool& operator=(const WavefrontPool&) = delete;
    
    // Allow move semantics (e.g., for transferring context between threads)
    WavefrontPool(WavefrontPool&&) noexcept = default;
    WavefrontPool& operator=(WavefrontPool&&) noexcept = default;

    // =========================================================
    // 1. Core API
    // =========================================================
    
    // Acquires a pointer to a WavefrontArray from the pool.
    // Specifying expected_capacity > 0 prevents internal vector reallocation upon acquisition.
    WavefrontArray* acquire(size_t expected_capacity = 0) {
        WavefrontArray* wf = nullptr;
        
        if (_cursor < _pool.size()) {
            wf = _pool[_cursor].get();
        } else {
            // Dynamically expand if the pool is exhausted
            _pool.push_back(std::make_unique<WavefrontArray>());
            wf = _pool.back().get();
        }

        // Always reset the logical size to 0 to hide previous garbage data
        wf->clear_logical_size();

        if (expected_capacity > 0) {
            wf->reserve_capacity(expected_capacity);
        }

        _cursor++;
        return wf;
    }

    // Called upon completion of a query or alignment process to return all arrays to the "free" state.
    // Internal memory (Capacity) is not freed and is reused in the next execution.
    void release_all() {
        for (size_t i = 0; i < _cursor; ++i) {
            _pool[i]->clear_logical_size();
        }
        _cursor = 0;
    }

    // =========================================================
    // 2. Diagnostics API
    // =========================================================
    
    // Total number of WavefrontArrays currently allocated in the pool
    size_t capacity() const {
        return _pool.size();
    }

    // Number of WavefrontArrays currently in use (loaned out)
    size_t used_count() const {
        return _cursor;
    }
    
    // Physically frees all allocated memory to resolve memory pressure
    void free_memory() {
        _pool.clear();
        _pool.shrink_to_fit();
        _cursor = 0;
    }
};
