#include "PatriciaWFA.hpp"

// =========================================================
// テンプレート関数の実装
// =========================================================
template <typename CostType>
template <typename StopPredicate>
std::vector<AlignmentResult> PatriciaWFA<CostType>::search_kernel(
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

    std::vector<bool> found_string_ids(_patricia_tree.size(), false);
    
    size_t history_size;
    if constexpr (CostType::is_unit) {  // Simple Edit Distances
        history_size = 2;
    } else {  // Linear
        history_size = std::max(_cost.mismatch, _cost.gap) + 1;
    }

    std::vector<WavefrontArray> wf_history(history_size);
    // WavefrontArray wf_array;
    WavefrontArray next_wf_array;
    WavefrontArray child_wf_array;
    std::array<WavefrontArray, PatriciaTree::CODE_MAX> buffer;
    
    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    uint32_t root = _patricia_tree.root_id();
    wf_history[0].push_back_state(root, 0, -1);
    
    int32_t current_score = 0;  // 現在の編集距離
    
    // Algorithm 1: GwfEditDist のメインループ
    while (true) {

        size_t curr_idx = current_score % history_size;
        WavefrontArray &curr_wf = wf_history[curr_idx];

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

        // Algorithm 2: GwfExtend
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
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(curr_vk);
            int32_t k = WavefrontArray::calc_k_from_vk(curr_vk);
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
                            if (found_string_ids[id]) {
                                continue;
                            }
                            found_string_ids[id] = true;
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

        // Algorithm 3: GwfExpand
        expand(padded_query, wf_history, next_wf_array, curr_idx, history_size, active_counts);

        if (upper_bound >= 0) {
            prune_by_upper_bound(
                wf_history[(current_score + 1) % history_size],
                subtree_max_lengths,
                subtree_min_lengths,
                query_length,
                upper_bound - (current_score + 1)
            );
        }

        ++current_score;

        // 以降のループで使うことはないので履歴をリセット
        wf_history[(current_score + 1) % history_size].clear_logical_size();
    }
    return results;
}

template <typename CostType>
template <typename StopPredicate>
std::vector<AlignmentResult> PatriciaWFA<CostType>::search_kernel(
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
    
    std::vector<bool> found_string_ids(_patricia_tree.size(), false);
    
    size_t history_size = std::max({ _cost.mismatch, _cost.gap_open + _cost.gap_extend, _cost.gap_extend }) + 1;

    std::vector<WavefrontArray> wf_history_d(history_size);
    std::vector<WavefrontArray> wf_history_m(history_size);
    std::vector<WavefrontArray> wf_history_i(history_size);
    // WavefrontArray wf_array;
    WavefrontArray next_wf_array_d;
    WavefrontArray next_wf_array_m;
    WavefrontArray next_wf_array_i;
    WavefrontArray child_wf_array;
    std::array<WavefrontArray, PatriciaTree::CODE_MAX> buffer;

    WavefrontArray pending_d;
    WavefrontArray merged_wf_array_d;
    
    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    uint32_t root = _patricia_tree.root_id();
    wf_history_m[0].push_back_state(root, 0, -1);
    
    int32_t current_score = 0;  // 現在の編集距離
    
    // Algorithm 1: GwfEditDist のメインループ
    while (true) {
        size_t curr_idx = current_score % history_size;
        WavefrontArray &curr_wf_m = wf_history_m[curr_idx];
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

        // Algorithm 2: GwfExtend
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
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(curr_vk);
            int32_t k = WavefrontArray::calc_k_from_vk(curr_vk);
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
                            if (found_string_ids[id]) {
                                continue;
                            }
                            found_string_ids[id] = true;
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

        // Algorithm 3: GwfExpand
        expand(padded_query, wf_history_d, wf_history_m, wf_history_i,
            next_wf_array_d, next_wf_array_m, next_wf_array_i,
            curr_idx, history_size, active_counts, buffer, pending_d, merged_wf_array_d);
        if (upper_bound >= 0) {
            prune_by_upper_bound<false>(
                wf_history_d[(current_score + 1) % history_size],
                wf_history_m[(current_score + 1) % history_size],
                wf_history_i[(current_score + 1) % history_size],
                _patricia_tree.get_subtree_max_lengths(),
                _patricia_tree.get_subtree_min_lengths(),
                query_length,
                upper_bound - (current_score + 1)
            );
        }

        ++current_score;

        // 以降のループで使うことはないので履歴をリセット
        wf_history_d[(current_score + 1) % history_size].clear_logical_size();
        wf_history_m[(current_score + 1) % history_size].clear_logical_size();
        wf_history_i[(current_score + 1) % history_size].clear_logical_size();
    }
    return results;
}
