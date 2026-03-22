#pragma once
#include <vector>
#include <cstdint>
#include <cstddef>

/**
 * @brief Avalanche mixer based on SplitMix64.
 * * A fast hash function capable of effectively scattering highly regular, 
 * consecutive integers (like vk keys) to prevent primary clustering 
 * in open-addressing hash tables.
 * * Reference:
 * Steele Jr, G. L., Lea, D., & Flood, C. H. (2014). 
 * "Fast splittable pseudorandom number generators." OOPSLA 2014.
 */
struct VkHash {
    inline size_t operator()(uint64_t vk) const noexcept {
        vk ^= vk >> 30;
        vk *= 0xbf58476d1ce4e5b9ULL;
        vk ^= vk >> 27;
        vk *= 0x94d049bb133111ebULL;
        vk ^= vk >> 31;
        return static_cast<size_t>(vk);
    }
};

template <typename Hash = VkHash>
class DynamicEpochHashMap {
private:
    struct Entry {
        uint64_t key;
        int32_t max_value;
        uint32_t epoch;
    };

    std::vector<Entry> _table;
    size_t _mask;
    uint32_t _current_epoch = 1;
    size_t _active_count = 0;
    
    Hash _hasher;

    void rehash() {
        size_t new_capacity = _table.size() * 2;
        // -2 represents an unreached/invalid offset.
        std::vector<Entry> new_table(new_capacity, {0, -2, 0});
        size_t new_mask = new_capacity - 1;

        for (const auto& e : _table) {
            if (e.epoch == _current_epoch) {
                size_t idx = _hasher(e.key) & new_mask;
                while (new_table[idx].epoch == _current_epoch) {
                    idx = (idx + 1) & new_mask;
                }
                new_table[idx] = e;
            }
        }
        _table = std::move(new_table);
        _mask = new_mask;
    }

public:
    explicit DynamicEpochHashMap(size_t initial_capacity_power_of_two = 65536) {
        _table.resize(initial_capacity_power_of_two, {0, -2, 0});
        _mask = initial_capacity_power_of_two - 1;
    }
    ~DynamicEpochHashMap() = default;

    DynamicEpochHashMap(const DynamicEpochHashMap&) = delete;
    DynamicEpochHashMap& operator=(const DynamicEpochHashMap&) = delete;

    DynamicEpochHashMap(DynamicEpochHashMap&&) noexcept = default;
    DynamicEpochHashMap& operator=(DynamicEpochHashMap&&) noexcept = default;

    inline bool update_and_check(uint64_t key, int32_t new_value) {
        // Rehash if load factor exceeds 5/8 (62.5%)
        if (_active_count > _table.size() * 5 / 8) {
            rehash();
        }

        size_t idx = _hasher(key) & _mask;

        while (true) {
            Entry& e = _table[idx];

            if (e.epoch != _current_epoch || e.key == key) {
                if (e.epoch == _current_epoch && new_value <= e.max_value) {
                    return false; 
                }
                
                if (e.epoch != _current_epoch) {
                    _active_count++;
                }
                
                e.key = key;
                e.max_value = new_value;
                e.epoch = _current_epoch;
                return true;
            }
            idx = (idx + 1) & _mask;
        }
    }

    void reset() {
        _current_epoch++;
        _active_count = 0;
    }
};
