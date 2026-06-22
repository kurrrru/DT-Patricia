#pragma once

#include "PatriciaTree.hpp"
#include "AlignmentResult.hpp"
#include "WavefrontArray.hpp"
#include "FastLCP.hpp"
#include "CostType.hpp"
#include "mod_utils.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <string_view>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <array>

template <typename CostType = UnitCost>
class DTPatricia {
 public:
    // =========================================================
    // 1. コンストラクタ / デストラクタ (Rule of Five)
    // =========================================================

    DTPatricia(const PatriciaTree &patricia_tree, CostType cost = CostType())
        : _patricia_tree(patricia_tree), _cost(cost) {}
    DTPatricia() = delete;
    ~DTPatricia() = default;
    DTPatricia(const DTPatricia&) = delete;
    DTPatricia &operator=(const DTPatricia&) = delete;
    DTPatricia(DTPatricia&&) noexcept = delete;
    DTPatricia &operator=(DTPatricia&&) noexcept = delete;

    // =========================================================
    // 2. 基本API
    // =========================================================
    
    [[nodiscard]] inline const PatriciaTree &get_patricia_tree() const noexcept {
        return _patricia_tree;
    }

    // =========================================================
    // 3. アラインメントAPI
    // =========================================================

    std::vector<AlignmentResult> ed_to_all(const std::string &query) const {
        std::size_t max_results = _patricia_tree.string_count();
        return search_kernel(query, [max_results](int, const auto& r){ return r.size() == max_results; }, -1);
    }

    std::vector<AlignmentResult> ed_within_k(const std::string &query, int k) const {
        std::size_t max_results = _patricia_tree.string_count();
        return search_kernel(query, [k, max_results](int current_score, const auto& r){ return current_score >= k || r.size() == max_results; }, k);
    }

    std::vector<AlignmentResult> ed_kth_smallest(const std::string &query, size_t k) const {
        if (k > _patricia_tree.string_count()) {
            k = _patricia_tree.string_count();
        }
        return search_kernel(query, [k](int, const auto& r){ return r.size() >= k; }, -1);
    }

    template <typename StopPredicate>
    std::vector<AlignmentResult> search_kernel(const std::string &query, StopPredicate stop_predicate, int upper_bound = -1) const requires (CostType::is_linear);

    template <typename StopPredicate>
    std::vector<AlignmentResult> search_kernel(const std::string &query, StopPredicate stop_predicate, int upper_bound = -1) const requires (!CostType::is_linear);

 private:
    const PatriciaTree &_patricia_tree;
    CostType _cost;

    void prune_by_upper_bound(WavefrontArray &wf_array,
                const std::vector<uint32_t> &subtree_max_lengths,
                const std::vector<uint32_t> &subtree_min_lengths,
                int32_t query_length,
                int upper_bound_remain
                ) const requires (CostType::is_linear);

    template <bool only_m>
    void prune_by_upper_bound(WavefrontArray &wf_array_d,
                WavefrontArray &wf_array_m,
                WavefrontArray &wf_array_i,
                const std::vector<uint32_t> &subtree_max_lengths,
                const std::vector<uint32_t> &subtree_min_lengths,
                int32_t query_length,
                int upper_bound_remain
                ) const requires (!CostType::is_linear);

    void extend(const std::string_view query, WavefrontArray &wf_array,
                    WavefrontArray &next_wf_array,
                    WavefrontArray &child_wf_array,
                    std::array<WavefrontArray, 5> &buffer,
                    const std::vector<uint32_t> &active_counts
                ) const;
    
    void expand(const std::string_view query,
                    std::vector<WavefrontArray> &wf_history,
                    WavefrontArray &next_wf_array,
                    int32_t curr_idx,
                    size_t history_size,
                    const std::vector<uint32_t> &active_counts,
                    std::vector<int32_t> &expand_scratch
                ) const requires (CostType::is_linear);

    void expand(const std::string_view query,
                    std::vector<WavefrontArray> &wf_history_d,
                    std::vector<WavefrontArray> &wf_history_m,
                    std::vector<WavefrontArray> &wf_history_i,
                    WavefrontArray &next_wf_array_d,
                    WavefrontArray &next_wf_array_m,
                    WavefrontArray &next_wf_array_i,
                    int32_t curr_idx,
                    size_t history_size,
                    const std::vector<uint32_t> &active_counts,
                    std::array<WavefrontArray, 5> &pending_d_buffer,
                    WavefrontArray &pending_d,
                    WavefrontArray &merged_wf_array_d,
                    std::vector<int32_t> &expand_scratch
                ) const requires (!CostType::is_linear);
};

#include "DT-Patricia_search_kernel.tpp"
#include "DT-Patricia_extend.tpp"
#include "DT-Patricia_expand.tpp"
#include "DT-Patricia_pruning.tpp"