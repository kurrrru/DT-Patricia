#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// =========================================================
// テンプレート関数の実装
// =========================================================
template <typename CostType>
template <typename StopPredicate>
std::vector<AlignmentResult> DTPatricia<CostType>::search_kernel(
    const std::string &query,
    StopPredicate stop_predicate,
    int upper_bound)
    const requires (CostType::is_linear) {

    std::vector<AlignmentResult> results;

    std::vector<uint32_t> active_counts = _patricia_tree.get_subtree_counts();
    const std::vector<uint32_t>& subtree_max_lengths = _patricia_tree.get_subtree_max_lengths();
    const std::vector<uint32_t>& subtree_min_lengths = _patricia_tree.get_subtree_min_lengths();

    if (_patricia_tree.empty() || query.empty()) {
        return results;
    }

    // =========================================================
    // クエリのパディング処理
    // =========================================================
    std::string padded_query_str = query;
    padded_query_str.append(PatriciaTree::SIMD_PADDING_SIZE, '\0');
    std::string_view padded_query(padded_query_str.data(), query.length());
    const int32_t query_length = static_cast<int32_t>(padded_query.length());

    std::vector<uint8_t> found_string_ids(_patricia_tree.string_count(), 0);

    uint32_t history_size;
    if constexpr (CostType::is_unit) {  // Simple Edit Distances
        history_size = 2;
    } else {  // Linear
        history_size = std::max(_cost.mismatch, _cost.gap) + 1;
    }

    std::vector<internal::WavefrontArray> wf_history(history_size);
    internal::WavefrontArray next_wf_array;
    internal::WavefrontArray child_wf_array;
    std::array<internal::WavefrontArray, PatriciaTree::CODE_MAX> buffer;
    std::vector<int32_t> expand_scratch;

    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    uint32_t root = _patricia_tree.root_id();
    wf_history[0].push_back_state(root, 0, -1);

    int32_t current_score = 0;  // 現在の編集距離
    uint32_t curr_idx = 0;

    // Algorithm 1: DT-Patricia のメインループ
    while (true) {
        internal::WavefrontArray &curr_wf = wf_history[curr_idx];

        if (curr_wf.empty()) {
            bool any_nonempty = false;
            for (size_t h = 0; h < history_size; ++h) {
                if (!wf_history[h].empty()) {
                    any_nonempty = true;
                    break;
                }
            }
            if (!any_nonempty) {
                break;  // すべての履歴が空 → 終了
            }
        }

        // Algorithm 2: DT-Patricia Extend
        extend(padded_query, curr_wf, next_wf_array, child_wf_array, buffer, active_counts);

        if (upper_bound >= 0) {
            prune_by_upper_bound(
                curr_wf,
                subtree_max_lengths,
                subtree_min_lengths,
                query_length,
                upper_bound - current_score
            );
        }

        // 終端チェック: クエリ全体が処理されたノードを探す
        for (size_t idx = 0; idx < curr_wf.active_size(); ++idx) {
            uint64_t curr_vk = curr_wf.get_vk(idx);
            uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vk(curr_vk);
            int32_t k = internal::WavefrontArray::calc_k_from_vk(curr_vk);
            int32_t j = curr_wf.get_offset(idx);
            int32_t i = k + j;

            if (i + 1 == query_length) {  // クエリ終端に到達
                if (j + 1 == static_cast<int32_t>(_patricia_tree.get_label_length(node_id))) {
                    // 終端文字列を持つかチェック
                    uint32_t term_node = _patricia_tree.transition(node_id, PatriciaTree::CODE_TERM);
                    if (term_node != 0) {
                        auto string_ids = _patricia_tree.get_string_id(term_node);
                        uint32_t found_count = 0;
                        for (uint32_t id : string_ids) {
                            if (found_string_ids[id] != 0) {
                                continue;
                            }
                            found_string_ids[id] = 1;
                            results.emplace_back(id, current_score);
                            found_count++;
                        }
                        if (found_count > 0) {
                            uint32_t curr_node = node_id;
                            while (curr_node != 0) {
                                if (active_counts[curr_node] > 0) {
                                    active_counts[curr_node] -= found_count;
                                }
                                curr_node = _patricia_tree.get_parent(curr_node);
                            }
                        }
                    }
                }
            }
        }

        // 停止条件チェック
        if (stop_predicate(current_score, results)) {
            break;
        }

        // Algorithm 3: DT-Patricia Expand
        expand(padded_query, wf_history, next_wf_array, curr_idx, history_size, active_counts, expand_scratch);

        uint32_t next_idx = internal::increment_mod(curr_idx, history_size);
        if (upper_bound >= 0) {
            prune_by_upper_bound(
                wf_history[next_idx],
                subtree_max_lengths,
                subtree_min_lengths,
                query_length,
                upper_bound - (current_score + 1)
            );
        }

        ++current_score;
        curr_idx = next_idx;
        next_idx = internal::increment_mod(next_idx, history_size);

        // 以降のループで使うことはないので履歴をリセット
        wf_history[next_idx].clear_logical_size();
    }
    return results;
}

template <typename CostType>
template <typename StopPredicate>
std::vector<AlignmentResult> DTPatricia<CostType>::search_kernel(
    const std::string &query,
    StopPredicate stop_predicate,
    int upper_bound)
    const requires (!CostType::is_linear) {
    (void)upper_bound;  // 未使用

    std::vector<AlignmentResult> results;

    std::vector<uint32_t> active_counts = _patricia_tree.get_subtree_counts();

    if (_patricia_tree.empty() || query.empty()) {
        return results;
    }

    // =========================================================
    // クエリのパディング処理
    // =========================================================
    std::string padded_query_str = query;
    padded_query_str.append(PatriciaTree::SIMD_PADDING_SIZE, '\0');
    std::string_view padded_query(padded_query_str.data(), query.length());
    const int32_t query_length = static_cast<int32_t>(padded_query.length());

    std::vector<uint8_t> found_string_ids(_patricia_tree.string_count(), 0);

    uint32_t history_size = std::max({ _cost.mismatch, _cost.gap_open + _cost.gap_extend, _cost.gap_extend }) + 1;

    std::vector<internal::WavefrontArray> wf_history_d(history_size);
    std::vector<internal::WavefrontArray> wf_history_m(history_size);
    std::vector<internal::WavefrontArray> wf_history_i(history_size);
    internal::WavefrontArray next_wf_array_d;
    internal::WavefrontArray next_wf_array_m;
    internal::WavefrontArray next_wf_array_i;
    internal::WavefrontArray child_wf_array;
    std::array<internal::WavefrontArray, PatriciaTree::CODE_MAX> buffer;

    internal::WavefrontArray pending_d;
    internal::WavefrontArray merged_wf_array_d;
    std::vector<int32_t> expand_scratch;

    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    uint32_t root = _patricia_tree.root_id();
    wf_history_m[0].push_back_state(root, 0, -1);

    int32_t current_score = 0;  // 現在の編集距離
    uint32_t curr_idx = 0;

    // Algorithm 1: DT-Patricia のメインループ
    while (true) {
        internal::WavefrontArray &curr_wf_m = wf_history_m[curr_idx];
        if (curr_wf_m.empty()) {
            bool any_nonempty = false;
            for (size_t h = 0; h < history_size; ++h) {
                if (!wf_history_m[h].empty()) {
                    any_nonempty = true;
                    break;
                }
            }
            if (!any_nonempty) {
                for (size_t h = 0; h < history_size; ++h) {
                    if (!wf_history_d[h].empty()) {
                        any_nonempty = true;
                        break;
                    }
                }
            }
            if (!any_nonempty) {
                for (size_t h = 0; h < history_size; ++h) {
                    if (!wf_history_i[h].empty()) {
                        any_nonempty = true;
                        break;
                    }
                }
            }
            if (!any_nonempty) {
                break;  // すべての履歴が空 → 終了
            }
        }

        // Algorithm 2: DT-Patricia Extend
        extend(padded_query, curr_wf_m, next_wf_array_m, child_wf_array, buffer, active_counts);
        if (upper_bound >= 0) {
            prune_by_upper_bound<true>(
                next_wf_array_d,
                curr_wf_m,
                next_wf_array_i,
                _patricia_tree.get_subtree_max_lengths(),
                _patricia_tree.get_subtree_min_lengths(),
                query_length,
                upper_bound - current_score
            );
        }

        // 終端チェック: クエリ全体が処理されたノードを探す
        for (size_t idx = 0; idx < curr_wf_m.active_size(); ++idx) {
            uint64_t curr_vk = curr_wf_m.get_vk(idx);
            uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vk(curr_vk);
            int32_t k = internal::WavefrontArray::calc_k_from_vk(curr_vk);
            int32_t j = curr_wf_m.get_offset(idx);
            int32_t i = k + j;

            if (i + 1 == query_length) {  // クエリ終端に到達
                if (j + 1 == static_cast<int32_t>(_patricia_tree.get_label_length(node_id))) {
                    // 終端文字列を持つかチェック
                    uint32_t term_node = _patricia_tree.transition(node_id, PatriciaTree::CODE_TERM);
                    if (term_node != 0) {
                        auto string_ids = _patricia_tree.get_string_id(term_node);
                        uint32_t found_count = 0;
                        for (uint32_t id : string_ids) {
                            if (found_string_ids[id] != 0) {
                                continue;
                            }
                            found_string_ids[id] = 1;
                            results.emplace_back(id, current_score);
                            found_count++;
                        }
                        if (found_count > 0) {
                            uint32_t curr_node = node_id;
                            while (curr_node != 0) {
                                if (active_counts[curr_node] > 0) {
                                    active_counts[curr_node] -= found_count;
                                }
                                curr_node = _patricia_tree.get_parent(curr_node);
                            }
                        }
                    }
                }
            }
        }

        // 停止条件チェック
        if (stop_predicate(current_score, results)) {
            break;
        }

        // Algorithm 3: DT-Patricia Expand
        expand(padded_query, wf_history_d, wf_history_m, wf_history_i,
            next_wf_array_d, next_wf_array_m, next_wf_array_i,
            curr_idx, history_size, active_counts, buffer, pending_d, merged_wf_array_d, expand_scratch);
        uint32_t next_idx = internal::increment_mod(curr_idx, history_size);
        if (upper_bound >= 0) {
            prune_by_upper_bound<false>(
                wf_history_d[next_idx],
                wf_history_m[next_idx],
                wf_history_i[next_idx],
                _patricia_tree.get_subtree_max_lengths(),
                _patricia_tree.get_subtree_min_lengths(),
                query_length,
                upper_bound - (current_score + 1)
            );
        }

        ++current_score;
        curr_idx = next_idx;
        next_idx = internal::increment_mod(next_idx, history_size);

        // 以降のループで使うことはないので履歴をリセット
        wf_history_d[next_idx].clear_logical_size();
        wf_history_m[next_idx].clear_logical_size();
        wf_history_i[next_idx].clear_logical_size();
    }
    return results;
}

}  // namespace dt_patricia
