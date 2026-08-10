#pragma once

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <queue>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <dt_patricia/alignment_result.hpp>
#include <dt_patricia/internal/lcp.hpp>
#include <dt_patricia/internal/mod_arithmetic.hpp>
#include <dt_patricia/internal/reached_offset_table.hpp>
#include <dt_patricia/internal/wavefront_array.hpp>
#include <dt_patricia/patricia_tree.hpp>
#include <dt_patricia/policy/cost.hpp>

namespace dt_patricia {

template <AlphabetPolicy Alphabet = DnaAlphabet, typename CostType = UnitCost>
class DTPatricia {
 public:
    using alphabet_type = Alphabet;
    using tree_type = PatriciaTree<Alphabet>;

    // =========================================================
    // 1. コンストラクタ / デストラクタ (Rule of Five)
    // =========================================================

    DTPatricia(const tree_type &patricia_tree, CostType cost = CostType())
        : _patricia_tree(patricia_tree), _cost(cost) {}
    DTPatricia() = delete;
    ~DTPatricia() = default;
    DTPatricia(const DTPatricia &) = delete;
    DTPatricia &operator=(const DTPatricia &) = delete;
    DTPatricia(DTPatricia &&) noexcept = delete;
    DTPatricia &operator=(DTPatricia &&) noexcept = delete;

    // =========================================================
    // 2. 基本API
    // =========================================================

    [[nodiscard]] inline const tree_type &get_patricia_tree() const noexcept {
        return _patricia_tree;
    }

    // =========================================================
    // 3. アラインメントAPI
    // =========================================================

    std::vector<AlignmentResult> ed_to_all(const std::string &query) const {
        std::size_t max_results = _patricia_tree.string_count();
        return search_kernel(
            query, [max_results](int, const auto &r) { return r.size() == max_results; }, -1);
    }

    std::vector<AlignmentResult> ed_within_k(const std::string &query, int k) const {
        std::size_t max_results = _patricia_tree.string_count();
        return search_kernel(
            query,
            [k, max_results](int current_score, const auto &r) {
                return current_score >= k || r.size() == max_results;
            },
            k);
    }

    std::vector<AlignmentResult> ed_kth_smallest(const std::string &query, size_t k) const {
        if (k > _patricia_tree.string_count()) {
            k = _patricia_tree.string_count();
        }
        return search_kernel(query, [k](int, const auto &r) { return r.size() >= k; }, -1);
    }

    template <typename StopPredicate>
    std::vector<AlignmentResult> search_kernel(const std::string &query,
                                               StopPredicate stop_predicate,
                                               int upper_bound = -1) const
        requires(CostType::is_linear);

    template <typename StopPredicate>
    std::vector<AlignmentResult> search_kernel(const std::string &query,
                                               StopPredicate stop_predicate,
                                               int upper_bound = -1) const
        requires(!CostType::is_linear);

 private:
    const tree_type &_patricia_tree;
    CostType _cost;

    void prune_by_upper_bound(internal::WavefrontArray &wf_array,
                              const std::vector<uint32_t> &subtree_max_lengths,
                              const std::vector<uint32_t> &subtree_min_lengths,
                              int32_t query_length, int upper_bound_remain) const
        requires(CostType::is_linear);

    template <bool only_m>
    void prune_by_upper_bound(internal::WavefrontArray &wf_array_d,
                              internal::WavefrontArray &wf_array_m,
                              internal::WavefrontArray &wf_array_i,
                              const std::vector<uint32_t> &subtree_max_lengths,
                              const std::vector<uint32_t> &subtree_min_lengths,
                              int32_t query_length, int upper_bound_remain) const
        requires(!CostType::is_linear);

    void extend(const std::string_view query, internal::WavefrontArray &wf_array,
                internal::WavefrontArray &next_wf_array, internal::WavefrontArray &child_wf_array,
                std::array<internal::WavefrontArray, PatriciaTree<Alphabet>::CODE_MAX> &buffer,
                const std::vector<uint32_t> &active_counts, internal::ReachedOffsetTable &reached,
                int32_t current_score) const;

    // expand_maxj は UnitCost の密経路で、出力対角線ごとの max_j を一旦受けるための作業領域。
    // 計算と書き出しを分けることで、計算側を分岐のない要素ごとのループにしている。
    void expand(const std::string_view query, std::vector<internal::WavefrontArray> &wf_history,
                internal::WavefrontArray &next_wf_array, int32_t curr_idx, size_t history_size,
                const std::vector<uint32_t> &active_counts, std::vector<int32_t> &expand_scratch,
                std::vector<int32_t> &expand_maxj) const
        requires(CostType::is_linear);

    void expand(
        const std::string_view query, std::vector<internal::WavefrontArray> &wf_history_d,
        std::vector<internal::WavefrontArray> &wf_history_m,
        std::vector<internal::WavefrontArray> &wf_history_i,
        internal::WavefrontArray &next_wf_array_d, internal::WavefrontArray &next_wf_array_m,
        internal::WavefrontArray &next_wf_array_i, int32_t curr_idx, size_t history_size,
        const std::vector<uint32_t> &active_counts,
        std::array<internal::WavefrontArray, PatriciaTree<Alphabet>::CODE_MAX> &pending_d_buffer,
        internal::WavefrontArray &pending_d, internal::WavefrontArray &merged_wf_array_d,
        std::vector<int32_t> &expand_scratch, internal::ReachedOffsetTable &reached_d,
        int32_t current_score) const
        requires(!CostType::is_linear);
};

}  // namespace dt_patricia

#include <dt_patricia/internal/detail_aligner/expand.tpp>
#include <dt_patricia/internal/detail_aligner/extend.tpp>
#include <dt_patricia/internal/detail_aligner/pruning.tpp>
#include <dt_patricia/internal/detail_aligner/search_kernel.tpp>
