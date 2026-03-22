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

    DynamicEpochHashMap<> visited_map;
    
    std::vector<AlignmentResult> results;

    std::vector<uint32_t> active_counts = _patricia_tree.get_subtree_counts();
    const std::vector<uint32_t>& subtree_max_lengths = _patricia_tree.get_subtree_max_lengths();
    const std::vector<uint32_t>& subtree_min_lengths = _patricia_tree.get_subtree_min_lengths();

    if (_patricia_tree.empty() || query.empty()) {
        return results;
    }

    const int32_t query_length = static_cast<int32_t>(query.length());
    
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
    std::array<WavefrontArray, 5> buffer;
    
    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    uint32_t root = _patricia_tree.root_id();
    // wf_array.push_back(root, 0, -1);
    wf_history[0].push_back(root, 0, -1);
    
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
        extend(query, curr_wf, next_wf_array, child_wf_array, buffer, active_counts, visited_map);
        
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
        for (size_t idx = 0; idx < curr_wf.size(); ++idx) {
            const WavefrontArray::DiagState &state = curr_wf[idx];
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(state.vk);
            int32_t k = WavefrontArray::calc_k_from_vk(state.vk);
            int32_t j = state.offset;
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
        expand(query, wf_history, next_wf_array, curr_idx, history_size, active_counts, visited_map);

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
        wf_history[(current_score + 1) % history_size].reset();
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

    DynamicEpochHashMap<> visited_map_d;
    DynamicEpochHashMap<> visited_map_m;
    DynamicEpochHashMap<> visited_map_i;

    std::vector<AlignmentResult> results;

    std::vector<uint32_t> active_counts = _patricia_tree.get_subtree_counts();

    if (_patricia_tree.empty() || query.empty()) {
        return results;
    }

    const int32_t query_length = static_cast<int32_t>(query.length());
    
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
    std::array<WavefrontArray, 5> buffer;
    
    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    uint32_t root = _patricia_tree.root_id();
    wf_history_m[0].push_back(root, 0, -1);
    
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
        extend(query, curr_wf_m, next_wf_array_m, child_wf_array, buffer, active_counts, visited_map_m);
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
        for (size_t idx = 0; idx < curr_wf_m.size(); ++idx) {
            const WavefrontArray::DiagState &state = curr_wf_m[idx];
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(state.vk);
            int32_t k = WavefrontArray::calc_k_from_vk(state.vk);
            int32_t j = state.offset;
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
        expand(query, wf_history_d, wf_history_m, wf_history_i,
            next_wf_array_d, next_wf_array_m, next_wf_array_i,
            curr_idx, history_size, active_counts,
            visited_map_d, visited_map_m, visited_map_i);
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
        wf_history_d[(current_score + 1) % history_size].reset();
        wf_history_m[(current_score + 1) % history_size].reset();
        wf_history_i[(current_score + 1) % history_size].reset();
    }
    return results;
}

// =========================================================
// Algorithm 2: GwfExtend - Exact match extension
// =========================================================
template <typename CostType>
void PatriciaWFA<CostType>::extend(
    const std::string_view query,
    WavefrontArray &wf_array,
    WavefrontArray &next_wf_array,
    WavefrontArray &child_wf_array,
    std::array<WavefrontArray, 5> &buffer,
    const std::vector<uint32_t> &active_counts,
    DynamicEpochHashMap<> &visited_map) const {
    next_wf_array.reset();
    child_wf_array.reset();
    for (auto &buf : buffer) {
        buf.reset();
    }
    bool buffer_used = false;
    
    const int32_t query_length = static_cast<int32_t>(query.length());
    
    size_t wf_array_idx = 0;
    size_t child_idx = 0;
    uint32_t last_node_id = UINT32_MAX;
    
    while (wf_array_idx < wf_array.size() || child_idx < child_wf_array.size() || buffer_used) {
        bool wf_array_idx_increment = false;
        bool child_idx_increment = false;

        // === ステップ1: 次に処理する状態を選択 ===
        WavefrontArray::DiagState current_state;
        if (child_idx >= child_wf_array.size() && wf_array_idx >= wf_array.size()) {
            // bufferのみ残っている場合
            for (const auto& buf : buffer) {
                for (size_t i = 0; i < buf.size(); ++i) {
                    const auto& state = buf[i];
                    child_wf_array.push_back(WavefrontArray::calc_node_id_from_vk(state.vk),
                                        WavefrontArray::calc_k_from_vk(state.vk),
                                        state.offset);
                                        
                }
            }
            for (auto &buf : buffer) {
                buf.reset();
            }
            buffer_used = false;
            continue;
        }

        if (child_idx >= child_wf_array.size()) {
            current_state = wf_array[wf_array_idx];
            wf_array_idx_increment = true;
        } else if (wf_array_idx >= wf_array.size()) {
            current_state = child_wf_array[child_idx];
            child_idx_increment = true;
        } else {
            const auto wf_state = wf_array[wf_array_idx];
            const auto child_state = child_wf_array[child_idx];
            
            if (wf_state.vk < child_state.vk) {
                current_state = wf_state;
                wf_array_idx_increment = true;
            } else if (wf_state.vk > child_state.vk) {
                current_state = child_state;
                child_idx_increment = true;
            } else {
                // vk が同じ → offset が大きい方を採用
                current_state = (wf_state.offset >= child_state.offset) ? wf_state : child_state;
                wf_array_idx_increment = true;
                child_idx_increment = true;
            }
        }
        
        // === ステップ2: node_idが変わったらbufferを処理 ===
        uint32_t node_id = WavefrontArray::calc_node_id_from_vk(current_state.vk);
        
        if (node_id != last_node_id && buffer_used) {
            // bufferはソート済み
            // (BFS順により child_wf_arrayの末尾 < bufferの最小値 が保証される)
            for (const auto& buf : buffer) {
                for (size_t i = 0; i < buf.size(); ++i) {
                    const auto& state = buf[i];
                    child_wf_array.push_back(WavefrontArray::calc_node_id_from_vk(state.vk),
                                        WavefrontArray::calc_k_from_vk(state.vk),
                                        state.offset);
                }
            }
            for (auto &buf : buffer) {
                buf.reset();
            }
            buffer_used = false;

            // これによってchild_wf_arrayに追加された値の方がwf_arrayの値より小さいかもしれない
            // child_idx_incrementが立っていない = child_idxの値はチェックされていない可能性がある
            if (!child_idx_increment) {
                const auto child_state = child_wf_array[child_idx];
                if (child_state.vk < current_state.vk) {
                    current_state = child_state;
                    child_idx_increment = true;
                    wf_array_idx_increment = false;
                } else if (child_state.vk == current_state.vk) {
                    if (child_state.offset > current_state.offset) {
                        current_state = child_state;
                    }
                    child_idx_increment = true;
                    wf_array_idx_increment = true;
                }
                node_id = WavefrontArray::calc_node_id_from_vk(current_state.vk);
            }
        }
        last_node_id = node_id;
        
        // === ステップ3: Extension 処理 ===
        int32_t k = WavefrontArray::calc_k_from_vk(current_state.vk);
        int32_t j = current_state.offset;
        int32_t i = k + j;
        
        std::string_view label = _patricia_tree.get_label(node_id);
        const int32_t label_len = static_cast<int32_t>(label.length());
        const int32_t max_lcp_len = std::min(query_length - (i + 1), label_len - (j + 1));
        // exact match extension
        int32_t lcp_len = static_cast<int32_t>(fast_lcp(query.data() + i + 1, label.data() + j + 1, max_lcp_len));
        i += lcp_len;
        j += lcp_len;
        
        // === ステップ4: 子ノードへの遷移 or next_wf_array への追加 ===
        if (j + 1 == label_len) {
            // ノード終端に到達 → 子ノードをbufferに追加
            for (uint8_t code = 1; code <= 5; ++code) {
                uint32_t child = _patricia_tree.transition(node_id, code);
                if (child != 0 && active_counts[child] > 0) {
                    int32_t new_k = (i + 1) - 0;
                    buffer[code - 1].push_back(child, new_k, -1);
                    buffer_used = true;
                }
            }
            
            if (_patricia_tree.is_terminal(node_id)) {
                next_wf_array.push_back(node_id, k, j);
            }
        } else {
            // ノード内で停止  ここで順序逆転起きる可能性あり
            next_wf_array.push_back(node_id, k, j);
        }
        if (wf_array_idx_increment) {
            ++wf_array_idx;
        }
        if (child_idx_increment) {
            ++child_idx;
        }
    }
    // ここではbufferは空になっている(空になっていないとwhileループが継続するため)
    wf_array.swap(next_wf_array);
}

// =========================================================
// Algorithm 3: GwfExpand - Edit operations (I/D/S)
// =======================================================
template <typename CostType>
void PatriciaWFA<CostType>::expand(
    const std::string_view query,
    std::vector<WavefrontArray> &wf_history,
    WavefrontArray &next_wf_array,
    int32_t curr_idx,
    size_t history_size,
    const std::vector<uint32_t> &active_counts,
    DynamicEpochHashMap<> &visited_map) const requires (CostType::is_linear) {
    next_wf_array.reset();
    const int32_t query_length = static_cast<int32_t>(query.length());
    // current_score + 1 に対応する WavefrontArray を構築する
    if constexpr (CostType::is_unit) {
        WavefrontArray &wf_array = wf_history[curr_idx];
        for (size_t idx = 0; idx < wf_array.size();) {
            // [start_idx, end_idx): node_id が同じ区間
            size_t start_idx = idx;
            size_t end_idx = idx + 1;

            const auto &first_state = wf_array[idx];
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(first_state.vk);

            // node_idが等しい区間を探す
            while (end_idx < wf_array.size()) {
                if (WavefrontArray::calc_node_id_from_vk(wf_array[end_idx].vk) != node_id) {
                    break;
                }
                ++end_idx;
            }

            if (active_counts[node_id] == 0) {
                idx = end_idx;
                continue;
            }

            const int32_t label_len = static_cast<int32_t>(_patricia_tree.get_label_length(node_id));

            // 3つのストリーム (Delete, Subst, Insert) のインデックス
            // それぞれ wf_array[start_idx...end_idx) を走査する
            size_t idx_d = start_idx; // Deletion source (k -> k-1)
            size_t idx_s = start_idx; // Substitution source (k -> k)
            size_t idx_i = start_idx; // Insertion source (k -> k+1)

            // 3つのストリームがすべて処理し終わるまでループ
            while (idx_d < end_idx || idx_s < end_idx || idx_i < end_idx) {
                
                // 1. 次の出力 k の候補を探す (3つのうち最小の k を見つける)
                //    k_d_target = k_current - 1
                //    k_s_target = k_current
                //    k_i_target = k_current + 1
                
                int32_t k_d_target = (idx_d < end_idx) ? WavefrontArray::calc_k_from_vk(wf_array[idx_d].vk) - 1 : INT32_MAX;
                int32_t k_s_target = (idx_s < end_idx) ? WavefrontArray::calc_k_from_vk(wf_array[idx_s].vk)     : INT32_MAX;
                int32_t k_i_target = (idx_i < end_idx) ? WavefrontArray::calc_k_from_vk(wf_array[idx_i].vk) + 1 : INT32_MAX;

                int32_t min_k = std::min({k_d_target, k_s_target, k_i_target});
                
                int32_t max_j = INT32_MIN; // 無効値

                // 2. min_k に該当するストリームを進め、max_j を更新する
                if (min_k != INT32_MAX) {
                    // --- Deletion (from k = min_k + 1) ---
                    if (k_d_target == min_k) {
                        const auto &st = wf_array[idx_d];
                        if (st.offset + 1 < label_len) {
                            int32_t new_j = st.offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_d++;
                    }

                    // --- Substitution (from k = min_k) ---
                    if (k_s_target == min_k) {
                        const auto &st = wf_array[idx_s];
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                        if (i_pos + 1 < query_length && st.offset + 1 < label_len) {
                            int32_t new_j = st.offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_s++;
                    }

                    // --- Insertion (from k = min_k - 1) ---
                    if (k_i_target == min_k) {
                        const auto &st = wf_array[idx_i];
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                        // 条件: i + 1 < q_len
                        if (i_pos + 1 < query_length) {
                            int32_t new_j = st.offset;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_i++;
                    }

                    // 3. 結果の登録 (有効な更新があった場合のみ)
                    if (max_j >= -1) {
                        next_wf_array.push_back(node_id, min_k, max_j);
                    }
                }
            }
            idx = end_idx;
        }
        // 次の探索セルとして追加
        wf_history[(curr_idx + 1) % history_size].swap(next_wf_array);
    } else {
        size_t wf_history_idx_d = (curr_idx + 1 + history_size - _cost.gap) % history_size;
        size_t wf_history_idx_s = (curr_idx + 1 + history_size - _cost.mismatch) % history_size;
        size_t wf_history_idx_i = (curr_idx + 1 + history_size - _cost.gap) % history_size;

        WavefrontArray &wf_array_d = wf_history[wf_history_idx_d];
        WavefrontArray &wf_array_s = wf_history[wf_history_idx_s];
        WavefrontArray &wf_array_i = wf_history[wf_history_idx_i];

        size_t start_idx_d = 0; // Deletion source (k -> k - gap) wf_history[curr_idx - _cost.gap]を見る
        size_t start_idx_s = 0; // Substitution source (k -> k) wf_history[curr_idx - cost.mismatch]を見る
        size_t start_idx_i = 0; // Insertion source (k -> k + gap) wf_history[curr_idx - _cost.gap]を見る
        // wf_history[curr_idx - _cost.gap][start_idx_d], wf_history[curr_idx - cost.mismatch][start_idx_s], wf_history[curr_idx - _cost.gap][start_idx_i]の中でnode_idが最小のものを見つける
        // そのnode_idとnode_idが等しい区間をそれぞれ見つける (end_idx_d, end_idx_s, end_idx_iを決定する)
        // 以降、そのnode_id に対して単位コスト版と同様の処理を行う
        while (start_idx_d < wf_array_d.size() || start_idx_s < wf_array_s.size() || start_idx_i < wf_array_i.size()) {
            // 次のnode_idを決定
            uint32_t node_id = UINT32_MAX;
            if (start_idx_d < wf_array_d.size()) {
                uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_d[start_idx_d].vk);
                node_id = std::min(node_id, nid);
            }
            if (start_idx_s < wf_array_s.size()) {
                uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_s[start_idx_s].vk);
                node_id = std::min(node_id, nid);
            }
            if (start_idx_i < wf_array_i.size()) {
                uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_i[start_idx_i].vk);
                node_id = std::min(node_id, nid);
            }

            // node_idが等しい区間を見つける
            size_t end_idx_d = start_idx_d;
            while (end_idx_d < wf_array_d.size() &&
                   WavefrontArray::calc_node_id_from_vk(wf_array_d[end_idx_d].vk) == node_id) {
                ++end_idx_d;
            }
            size_t end_idx_s = start_idx_s;
            while (end_idx_s < wf_array_s.size() &&
                   WavefrontArray::calc_node_id_from_vk(wf_array_s[end_idx_s].vk) == node_id) {
                ++end_idx_s;
            }
            size_t end_idx_i = start_idx_i;
            while (end_idx_i < wf_array_i.size() &&
                   WavefrontArray::calc_node_id_from_vk(wf_array_i[end_idx_i].vk) == node_id) {
                ++end_idx_i;
            }

            if (active_counts[node_id] == 0) {
                start_idx_d = end_idx_d;
                start_idx_s = end_idx_s;
                start_idx_i = end_idx_i;
                continue;
            }

            const int32_t label_len = static_cast<int32_t>(_patricia_tree.get_label_length(node_id));

            // 3つのストリーム (Delete, Subst, Insert) のインデックス
            // それぞれ wf_array[start_idx...end_idx) を走査する
            size_t idx_d = start_idx_d; // Deletion source (k -> k-1)
            size_t idx_s = start_idx_s; // Substitution source (k -> k)
            size_t idx_i = start_idx_i; // Insertion source (k -> k+1)
            // 3つのストリームがすべて処理し終わるまでループ
            while (idx_d < end_idx_d || idx_s < end_idx_s || idx_i < end_idx_i) {
                // 1. 次の出力 k の候補を探す (3つのうち最小の k を見つける)
                //    k_d_target = k_current - 1
                //    k_s_target = k_current
                //    k_i_target = k_current + 1
                int32_t k_d_target = (idx_d < end_idx_d) ? WavefrontArray::calc_k_from_vk(wf_array_d[idx_d].vk) - 1 : INT32_MAX;
                int32_t k_s_target = (idx_s < end_idx_s) ? WavefrontArray::calc_k_from_vk(wf_array_s[idx_s].vk)     : INT32_MAX;
                int32_t k_i_target = (idx_i < end_idx_i) ? WavefrontArray::calc_k_from_vk(wf_array_i[idx_i].vk) + 1 : INT32_MAX;

                int32_t min_k = std::min({k_d_target, k_s_target, k_i_target});

                int32_t max_j = INT32_MIN; // 無効値

                if (min_k != INT32_MAX) {
                    // --- Deletion (from k = min_k + 1) ---
                    if (k_d_target == min_k) {
                        const auto &st = wf_array_d[idx_d];
                        if (st.offset + 1 < label_len) {
                            int32_t new_j = st.offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_d++;
                    }

                    // --- Substitution (from k = min_k) ---
                    if (k_s_target == min_k) {
                        const auto &st = wf_array_s[idx_s];
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                        if (i_pos + 1 < query_length && st.offset + 1 < label_len) {
                            int32_t new_j = st.offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_s++;
                    }

                    // --- Insertion (from k = min_k - 1) ---
                    if (k_i_target == min_k) {
                        const auto &st = wf_array_i[idx_i];
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                        // 条件: i + 1 < q_len
                        if (i_pos + 1 < query_length) {
                            int32_t new_j = st.offset;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_i++;
                    }

                    // 3. 結果の登録 (有効な更新があった場合のみ)
                    if (max_j >= -1) {
                        next_wf_array.push_back(node_id, min_k, max_j);
                    }
                }
            }
            start_idx_d = end_idx_d;
            start_idx_s = end_idx_s;
            start_idx_i = end_idx_i;
        }
        wf_history[(curr_idx + 1) % history_size].swap(next_wf_array);
    }
}

template <typename CostType>
void PatriciaWFA<CostType>::expand(
    const std::string_view query,
    std::vector<WavefrontArray> &wf_history_d,
    std::vector<WavefrontArray> &wf_history_m,
    std::vector<WavefrontArray> &wf_history_i,
    WavefrontArray &next_wf_array_d,
    WavefrontArray &next_wf_array_m,
    WavefrontArray &next_wf_array_i,
    int32_t curr_idx,
    size_t history_size, 
    const std::vector<uint32_t> &active_counts,
    DynamicEpochHashMap<> &visited_map_d,
    DynamicEpochHashMap<> &visited_map_m,
    DynamicEpochHashMap<> &visited_map_i) const requires (!CostType::is_linear) {
    next_wf_array_d.reset();
    next_wf_array_m.reset();
    next_wf_array_i.reset();
    const int32_t query_length = static_cast<int32_t>(query.length());
    
    // current_score + 1 に対応する WavefrontArray を構築する
    size_t wf_history_idx_d = (curr_idx + 1 + history_size - _cost.gap_extend) % history_size;
    size_t wf_history_idx_m = (curr_idx + 1 + history_size - _cost.gap_open - _cost.gap_extend) % history_size;
    size_t wf_history_idx_i = (curr_idx + 1 + history_size - _cost.gap_extend) % history_size;

    WavefrontArray &wf_array_d = wf_history_d[wf_history_idx_d];
    WavefrontArray &wf_array_m = wf_history_m[wf_history_idx_m];
    WavefrontArray &wf_array_i = wf_history_i[wf_history_idx_i];

    size_t start_idx_d = 0; 
    size_t start_idx_m = 0; 
    size_t start_idx_i = 0; 

    // ノード境界を越えるDeletion状態を一時保存するバッファ
    WavefrontArray pending_d;

    while (start_idx_d < wf_array_d.size() || start_idx_m < wf_array_m.size() || start_idx_i < wf_array_i.size()) {
        // 次のnode_idを決定
        uint32_t node_id = UINT32_MAX;
        if (start_idx_d < wf_array_d.size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_d[start_idx_d].vk);
            node_id = std::min(node_id, nid);
        }
        if (start_idx_m < wf_array_m.size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_m[start_idx_m].vk);
            node_id = std::min(node_id, nid);
        }
        if (start_idx_i < wf_array_i.size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_i[start_idx_i].vk);
            node_id = std::min(node_id, nid);
        }

        // node_idが等しい区間を見つける
        size_t end_idx_d = start_idx_d;
        while (end_idx_d < wf_array_d.size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_d[end_idx_d].vk) == node_id) {
            ++end_idx_d;
        }
        size_t end_idx_m = start_idx_m;
        while (end_idx_m < wf_array_m.size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_m[end_idx_m].vk) == node_id) {
            ++end_idx_m;
        }
        size_t end_idx_i = start_idx_i;
        while (end_idx_i < wf_array_i.size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_i[end_idx_i].vk) == node_id) {
            ++end_idx_i;
        }

        if (active_counts[node_id] == 0) {
            start_idx_d = end_idx_d;
            start_idx_m = end_idx_m;
            start_idx_i = end_idx_i;
            continue;
        }

        const int32_t label_len = static_cast<int32_t>(_patricia_tree.get_label_length(node_id));

        // 3つのストリーム (Delete, Subst, Insert) のインデックス
        size_t idx_d = start_idx_d; 
        size_t idx_m = start_idx_m; 
        
        while (idx_d < end_idx_d || idx_m < end_idx_m) {
            int32_t k_d_target = (idx_d < end_idx_d) ? WavefrontArray::calc_k_from_vk(wf_array_d[idx_d].vk) - 1 : INT32_MAX;
            int32_t k_m_target = (idx_m < end_idx_m) ? WavefrontArray::calc_k_from_vk(wf_array_m[idx_m].vk) - 1 : INT32_MAX;

            int32_t min_k = std::min({k_d_target, k_m_target});
            int32_t max_j = INT32_MIN; // 無効値

            if (min_k != INT32_MAX) {
                // --- Deletion (D -> D) ---
                if (k_d_target == min_k) {
                    const auto &st = wf_array_d[idx_d];
                    // ケース1: ノード内での伸長
                    if (st.offset + 1 < label_len) {
                        int32_t new_j = st.offset + 1;
                        max_j = std::max(max_j, new_j);
                    } 
                    // ケース2: ノード境界での遷移 (追加)
                    else {
                        int32_t current_k = WavefrontArray::calc_k_from_vk(st.vk);
                        // 次の k = i - next_j = (current_k + offset) - (-1) = current_k + offset + 1
                        int32_t next_k = current_k + st.offset;

                        for (uint8_t code = 1; code <= 5; ++code) {
                            uint32_t child = _patricia_tree.transition(node_id, code);
                            if (child != 0 && active_counts[child] > 0) {
                                pending_d.push_back(child, next_k, 0);
                            }
                        }
                    }
                    idx_d++;
                }

                // --- Deletion (M -> D) ---
                if (k_m_target == min_k) {
                    const auto &st = wf_array_m[idx_m];
                    if (st.offset + 1 < label_len) {
                        int32_t new_j = st.offset + 1;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_m++;
                }

                // 3. 結果の登録
                if (max_j >= -1) {
                    next_wf_array_d.push_back(node_id, min_k, max_j);
                }
            }
        }

        idx_m = start_idx_m;  
        size_t idx_i = start_idx_i;  
        
        while (idx_m < end_idx_m || idx_i < end_idx_i) {
            int32_t k_m_target = (idx_m < end_idx_m) ? WavefrontArray::calc_k_from_vk(wf_array_m[idx_m].vk) + 1 : INT32_MAX;
            int32_t k_i_target = (idx_i < end_idx_i) ? WavefrontArray::calc_k_from_vk(wf_array_i[idx_i].vk) + 1 : INT32_MAX;

            int32_t min_k = std::min({k_m_target, k_i_target});
            int32_t max_j = INT32_MIN; 

            if (min_k != INT32_MAX) {
                // --- Insertion (M -> I) ---
                if (k_m_target == min_k) {
                    const auto &st = wf_array_m[idx_m];
                    int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                    if (i_pos + 1 < query_length) {
                        int32_t new_j = st.offset;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_m++;
                }

                // --- Insertion (I -> I) ---
                if (k_i_target == min_k) {
                    const auto &st = wf_array_i[idx_i];
                    int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                    if (i_pos + 1 < query_length) {
                        int32_t new_j = st.offset;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_i++;
                }

                // 3. 結果の登録
                if (max_j >= -1) {
                    next_wf_array_i.push_back(node_id, min_k, max_j);
                }
            }
        }

        start_idx_d = end_idx_d;
        start_idx_m = end_idx_m;
        start_idx_i = end_idx_i;
    }

    // === 保留していたノード境界Deletionのマージ処理 ===
    if (!pending_d.empty()) {
        pending_d.dedup();

        WavefrontArray tmp_wf; 

        size_t idx1 = 0;
        size_t idx2 = 0;
        size_t len1 = next_wf_array_d.size();
        size_t len2 = pending_d.size();

        while (idx1 < len1 || idx2 < len2) {
            if (idx1 == len1) {
                const auto& st = pending_d[idx2];
                tmp_wf.push_back(WavefrontArray::calc_node_id_from_vk(st.vk), 
                                    WavefrontArray::calc_k_from_vk(st.vk), 
                                    st.offset);
                ++idx2;
            } else if (idx2 == len2) {
                const auto& st = next_wf_array_d[idx1];
                tmp_wf.push_back(WavefrontArray::calc_node_id_from_vk(st.vk), 
                                    WavefrontArray::calc_k_from_vk(st.vk), 
                                    st.offset);
                ++idx1;
            } else {
                const auto& st1 = next_wf_array_d[idx1];
                const auto& st2 = pending_d[idx2];

                if (st1.vk < st2.vk) {
                    tmp_wf.push_back(WavefrontArray::calc_node_id_from_vk(st1.vk), 
                                     WavefrontArray::calc_k_from_vk(st1.vk), 
                                     st1.offset);
                    ++idx1;
                } else if (st2.vk < st1.vk) {
                    tmp_wf.push_back(WavefrontArray::calc_node_id_from_vk(st2.vk), 
                                     WavefrontArray::calc_k_from_vk(st2.vk), 
                                     st2.offset);
                    ++idx2;
                } else {
                    const auto& winner = (st1.offset >= st2.offset) ? st1 : st2;
                    tmp_wf.push_back(WavefrontArray::calc_node_id_from_vk(winner.vk), 
                                        WavefrontArray::calc_k_from_vk(winner.vk), 
                                        winner.offset);
                    ++idx1;
                    ++idx2;
                }
            }
        }
        next_wf_array_d.swap(tmp_wf);
    }

    wf_history_d[(curr_idx + 1) % history_size].swap(next_wf_array_d);
    wf_history_i[(curr_idx + 1) % history_size].swap(next_wf_array_i);

    // match/mismatch のみを見る
    wf_history_idx_d = (curr_idx + 1) % history_size;
    wf_history_idx_m = (curr_idx + 1 + history_size - _cost.mismatch) % history_size;
    wf_history_idx_i = (curr_idx + 1) % history_size;

    WavefrontArray &wf_array_dm = wf_history_d[wf_history_idx_d];
    WavefrontArray &wf_array_mm = wf_history_m[wf_history_idx_m];
    WavefrontArray &wf_array_im = wf_history_i[wf_history_idx_i];

    start_idx_d = 0; 
    start_idx_m = 0; 
    start_idx_i = 0; 

    while (start_idx_d < wf_array_dm.size() || start_idx_m < wf_array_mm.size() || start_idx_i < wf_array_im.size()) {
        // 次のnode_idを決定
        uint32_t node_id = UINT32_MAX;
        if (start_idx_d < wf_array_dm.size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_dm[start_idx_d].vk);
            node_id = std::min(node_id, nid);
        }
        if (start_idx_m < wf_array_mm.size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_mm[start_idx_m].vk);
            node_id = std::min(node_id, nid);
        }
        if (start_idx_i < wf_array_im.size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_im[start_idx_i].vk);
            node_id = std::min(node_id, nid);
        }

        size_t end_idx_d = start_idx_d;
        while (end_idx_d < wf_array_dm.size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_dm[end_idx_d].vk) == node_id) {
            ++end_idx_d;
        }
        size_t end_idx_m = start_idx_m;
        while (end_idx_m < wf_array_mm.size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_mm[end_idx_m].vk) == node_id) {
            ++end_idx_m;
        }
        size_t end_idx_i = start_idx_i;
        while (end_idx_i < wf_array_im.size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_im[end_idx_i].vk) == node_id) {
            ++end_idx_i;
        }

        if (active_counts[node_id] == 0) {
            start_idx_d = end_idx_d;
            start_idx_m = end_idx_m;
            start_idx_i = end_idx_i;
            continue;
        }

        const int32_t label_len = static_cast<int32_t>(_patricia_tree.get_label_length(node_id));

        size_t idx_d = start_idx_d; 
        size_t idx_m = start_idx_m; 
        size_t idx_i = start_idx_i; 

        while (idx_d < end_idx_d || idx_m < end_idx_m || idx_i < end_idx_i) {
            int32_t k_d_target = (idx_d < end_idx_d) ? WavefrontArray::calc_k_from_vk(wf_array_dm[idx_d].vk)     : INT32_MAX;
            int32_t k_m_target = (idx_m < end_idx_m) ? WavefrontArray::calc_k_from_vk(wf_array_mm[idx_m].vk)     : INT32_MAX;
            int32_t k_i_target = (idx_i < end_idx_i) ? WavefrontArray::calc_k_from_vk(wf_array_im[idx_i].vk)     : INT32_MAX;

            int32_t min_k = std::min({k_d_target, k_m_target, k_i_target});
            int32_t max_j = INT32_MIN; 

            if (min_k != INT32_MAX) {
                if (k_d_target == min_k) {
                    const auto &st = wf_array_dm[idx_d];
                    int32_t new_j = st.offset;
                    max_j = std::max(max_j, new_j);
                    idx_d++;
                }

                if (k_m_target == min_k) {
                    const auto &st = wf_array_mm[idx_m];
                    int32_t i_pos = WavefrontArray::calc_k_from_vk(st.vk) + st.offset;
                    if (i_pos + 1 < query_length && st.offset + 1 < label_len) {
                        int32_t new_j = st.offset + 1;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_m++;
                }

                if (k_i_target == min_k) {
                    const auto &st = wf_array_im[idx_i];
                    int32_t new_j = st.offset;
                    max_j = std::max(max_j, new_j);
                    idx_i++;
                }

                if (max_j >= -1) {
                    next_wf_array_m.push_back(node_id, min_k, max_j);
                }
            }
        }
        start_idx_d = end_idx_d;
        start_idx_m = end_idx_m;
        start_idx_i = end_idx_i;
    }
    wf_history_m[(curr_idx + 1) % history_size].swap(next_wf_array_m);
}

// =========================================================
// prune_by_upper_bound
// =========================================================
template <typename CostType>
void PatriciaWFA<CostType>::prune_by_upper_bound(
    WavefrontArray &wf_array,
    const std::vector<uint32_t> &subtree_max_lengths,
    const std::vector<uint32_t> &subtree_min_lengths,
    int32_t query_length,
    int upper_bound_remain
    ) const requires (CostType::is_linear) {
    size_t write_idx = 0;
    for (size_t i = 0; i < wf_array.size(); i++) {
        const auto &state = wf_array[i];
        uint32_t node_id = WavefrontArray::calc_node_id_from_vk(state.vk);
        int32_t k = WavefrontArray::calc_k_from_vk(state.vk);
        int32_t j_pos = state.offset;
        int32_t i_pos = k + j_pos;
        int32_t max_remain = subtree_max_lengths[node_id] - (j_pos + 1);
        int32_t min_remain = subtree_min_lengths[node_id] - (j_pos + 1);
        int32_t query_remain = query_length - (i_pos + 1);
        if ((query_remain - max_remain) * static_cast<int32_t>(_cost.gap) <= upper_bound_remain &&
            (min_remain - query_remain) * static_cast<int32_t>(_cost.gap) <= upper_bound_remain) {
            if (write_idx != i) {
                wf_array[write_idx] = state;
            }
            write_idx++;
        }
    }
    wf_array.set_size(write_idx);
}

template <typename CostType>
template <bool only_m>
void PatriciaWFA<CostType>::prune_by_upper_bound(
    WavefrontArray &wf_array_d,
    WavefrontArray &wf_array_m,
    WavefrontArray &wf_array_i,
    const std::vector<uint32_t> &subtree_max_lengths,
    const std::vector<uint32_t> &subtree_min_lengths,
    int32_t query_length,
    int upper_bound_remain
) const requires (!CostType::is_linear) {

    const int32_t gap_o = static_cast<int32_t>(_cost.gap_open);
    const int32_t gap_e = static_cast<int32_t>(_cost.gap_extend);

    // M, I, D それぞれに対して判定を行う
    auto prune_logic = [&](WavefrontArray& wf, int state_type) {
        size_t write_idx = 0;
        for (size_t i = 0; i < wf.size(); ++i) {
            const auto &state = wf[i];
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(state.vk);
            int32_t v_k = WavefrontArray::calc_k_from_vk(state.vk);
            int32_t j_pos = state.offset;
            int32_t i_pos = v_k + j_pos;

            int32_t max_rem_t = static_cast<int32_t>(subtree_max_lengths[node_id]) - j_pos;
            int32_t min_rem_t = static_cast<int32_t>(subtree_min_lengths[node_id]) - j_pos;
            int32_t rem_q = query_length - i_pos;

            int32_t lb = 0;
            // クエリがターゲットより長い (Insertionが必要)
            if (rem_q > max_rem_t) {
                int32_t diff = rem_q - max_rem_t;
                // 現在が I 以外なら新たに gap_open が必要
                lb = (state_type == 2/*STATE_I*/) ? (diff * gap_e) : (gap_o + diff * gap_e);
            } 
            // ターゲットがクエリより長い (Deletionが必要)
            else if (min_rem_t > rem_q) {
                int32_t diff = min_rem_t - rem_q;
                // 現在が D 以外なら新たに gap_open が必要
                lb = (state_type == 1/*STATE_D*/) ? (diff * gap_e) : (gap_o + diff * gap_e);
            }

            if (lb <= upper_bound_remain) {
                if (write_idx != i) {
                    wf[write_idx] = state;
                }
                write_idx++;
            }
        }
        wf.set_size(write_idx);
    };

    prune_logic(wf_array_m, 0); // Match
    if constexpr (!only_m) {
        prune_logic(wf_array_d, 1); // Deletion
        prune_logic(wf_array_i, 2); // Insertion
    }
}

