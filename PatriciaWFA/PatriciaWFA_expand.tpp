
#include "PatriciaWFA.hpp"

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
    const std::vector<uint32_t> &active_counts) const requires (CostType::is_linear) {
    next_wf_array.clear_logical_size();
    const int32_t query_length = static_cast<int32_t>(query.length());
    // current_score + 1 に対応する WavefrontArray を構築する
    if constexpr (CostType::is_unit) {
        WavefrontArray &wf_array = wf_history[curr_idx];
        for (size_t idx = 0; idx < wf_array.active_size();) {
            // [start_idx, end_idx): node_id が同じ区間
            size_t start_idx = idx;
            size_t end_idx = idx + 1;

            const uint64_t first_vk = wf_array.get_vk(idx);
            uint32_t node_id = WavefrontArray::calc_node_id_from_vk(first_vk);

            // node_idが等しい区間を探す
            while (end_idx < wf_array.active_size()) {
                if (WavefrontArray::calc_node_id_from_vk(wf_array.get_vk(end_idx)) != node_id) {
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
                
                int32_t k_d_target = (idx_d < end_idx) ? WavefrontArray::calc_k_from_vk(wf_array.get_vk(idx_d)) - 1 : INT32_MAX;
                int32_t k_s_target = (idx_s < end_idx) ? WavefrontArray::calc_k_from_vk(wf_array.get_vk(idx_s))     : INT32_MAX;
                int32_t k_i_target = (idx_i < end_idx) ? WavefrontArray::calc_k_from_vk(wf_array.get_vk(idx_i)) + 1 : INT32_MAX;

                int32_t min_k = std::min({k_d_target, k_s_target, k_i_target});
                
                int32_t max_j = INT32_MIN; // 無効値

                // 2. min_k に該当するストリームを進め、max_j を更新する
                if (min_k != INT32_MAX) {
                    // --- Deletion (from k = min_k + 1) ---
                    if (k_d_target == min_k) {
                        const int32_t st_offset = wf_array.get_offset(idx_d);
                        if (st_offset + 1 < label_len) {
                            int32_t new_j = st_offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_d++;
                    }
                    
                    // --- Substitution (from k = min_k) ---
                    if (k_s_target == min_k) {
                        const int32_t st_offset = wf_array.get_offset(idx_s);
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(wf_array.get_vk(idx_s)) + st_offset;
                        if (i_pos + 1 < query_length && st_offset + 1 < label_len) {
                            int32_t new_j = st_offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_s++;
                    }

                    // --- Insertion (from k = min_k - 1) ---
                    if (k_i_target == min_k) {
                        const int32_t st_offset = wf_array.get_offset(idx_i);
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(wf_array.get_vk(idx_i)) + st_offset;
                        // 条件: i + 1 < q_len
                        if (i_pos + 1 < query_length) {
                            int32_t new_j = st_offset;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_i++;
                    }

                    // 3. 結果の登録 (有効な更新があった場合のみ)
                    if (max_j >= -1) {
                        next_wf_array.push_back_state(node_id, min_k, max_j);
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
        while (start_idx_d < wf_array_d.active_size() || start_idx_s < wf_array_s.active_size() || start_idx_i < wf_array_i.active_size()) {
            // 次のnode_idを決定
            uint32_t node_id = UINT32_MAX;
            if (start_idx_d < wf_array_d.active_size()) {
                uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(start_idx_d));
                node_id = std::min(node_id, nid);
            }
            if (start_idx_s < wf_array_s.active_size()) {
                uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_s.get_vk(start_idx_s));
                node_id = std::min(node_id, nid);
            }
            if (start_idx_i < wf_array_i.active_size()) {
                uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(start_idx_i));
                node_id = std::min(node_id, nid);
            }

            // node_idが等しい区間を見つける
            size_t end_idx_d = start_idx_d;
            while (end_idx_d < wf_array_d.active_size() &&
                   WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(end_idx_d)) == node_id) {
                ++end_idx_d;
            }
            size_t end_idx_s = start_idx_s;
            while (end_idx_s < wf_array_s.active_size() &&
                   WavefrontArray::calc_node_id_from_vk(wf_array_s.get_vk(end_idx_s)) == node_id) {
                ++end_idx_s;
            }
            size_t end_idx_i = start_idx_i;
            while (end_idx_i < wf_array_i.active_size() &&
                   WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(end_idx_i)) == node_id) {
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
                int32_t k_d_target = (idx_d < end_idx_d) ? WavefrontArray::calc_k_from_vk(wf_array_d.get_vk(idx_d)) - 1 : INT32_MAX;
                int32_t k_s_target = (idx_s < end_idx_s) ? WavefrontArray::calc_k_from_vk(wf_array_s.get_vk(idx_s))     : INT32_MAX;
                int32_t k_i_target = (idx_i < end_idx_i) ? WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) + 1 : INT32_MAX;

                int32_t min_k = std::min({k_d_target, k_s_target, k_i_target});

                int32_t max_j = INT32_MIN; // 無効値

                if (min_k != INT32_MAX) {
                    // --- Deletion (from k = min_k + 1) ---
                    if (k_d_target == min_k) {
                        // const auto &st = wf_array_d[idx_d];
                        const int32_t st_offset = wf_array_d.get_offset(idx_d);
                        if (st_offset + 1 < label_len) {
                            int32_t new_j = st_offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_d++;
                    }

                    // --- Substitution (from k = min_k) ---
                    if (k_s_target == min_k) {
                        const int32_t st_offset = wf_array_s.get_offset(idx_s);
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(wf_array_s.get_vk(idx_s)) + st_offset;
                        if (i_pos + 1 < query_length && st_offset + 1 < label_len) {
                            int32_t new_j = st_offset + 1;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_s++;
                    }

                    // --- Insertion (from k = min_k - 1) ---
                    if (k_i_target == min_k) {
                        const int32_t st_offset = wf_array_i.get_offset(idx_i);
                        int32_t i_pos = WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) + st_offset;
                        // 条件: i + 1 < q_len
                        if (i_pos + 1 < query_length) {
                            int32_t new_j = st_offset;
                            max_j = std::max(max_j, new_j);
                        }
                        idx_i++;
                    }

                    // 3. 結果の登録 (有効な更新があった場合のみ)
                    if (max_j >= -1) {
                        next_wf_array.push_back_state(node_id, min_k, max_j);
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
    std::array<WavefrontArray, PatriciaTree::CODE_MAX> &pending_d_buffer,
    WavefrontArray &pending_d,
    WavefrontArray &merged_wf_array_d) const requires (!CostType::is_linear) {
    next_wf_array_d.clear_logical_size();
    next_wf_array_m.clear_logical_size();
    next_wf_array_i.clear_logical_size();
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

    bool has_pending_d = false;

    while (start_idx_d < wf_array_d.active_size() || start_idx_m < wf_array_m.active_size() || start_idx_i < wf_array_i.active_size()) {
        // 次のnode_idを決定
        uint32_t node_id = UINT32_MAX;
        if (start_idx_d < wf_array_d.active_size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(start_idx_d));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_m < wf_array_m.active_size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_m.get_vk(start_idx_m));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_i < wf_array_i.active_size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(start_idx_i));
            node_id = std::min(node_id, nid);
        }

        // node_idが等しい区間を見つける
        size_t end_idx_d = start_idx_d;
        while (end_idx_d < wf_array_d.active_size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(end_idx_d)) == node_id) {
            ++end_idx_d;
        }
        size_t end_idx_m = start_idx_m;
        while (end_idx_m < wf_array_m.active_size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_m.get_vk(end_idx_m)) == node_id) {
            ++end_idx_m;
        }
        size_t end_idx_i = start_idx_i;
        while (end_idx_i < wf_array_i.active_size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(end_idx_i)) == node_id) {
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
            int32_t k_d_target = (idx_d < end_idx_d) ? WavefrontArray::calc_k_from_vk(wf_array_d.get_vk(idx_d)) - 1 : INT32_MAX;
            int32_t k_m_target = (idx_m < end_idx_m) ? WavefrontArray::calc_k_from_vk(wf_array_m.get_vk(idx_m)) - 1 : INT32_MAX;

            int32_t min_k = std::min({k_d_target, k_m_target});
            int32_t max_j = INT32_MIN; // 無効値

            if (min_k != INT32_MAX) {
                // --- Deletion (D -> D) ---
                if (k_d_target == min_k) {
                    const uint64_t st_vk = wf_array_d.get_vk(idx_d);
                    const int32_t st_offset = wf_array_d.get_offset(idx_d);
                    // ケース1: ノード内での伸長
                    if (st_offset + 1 < label_len) {
                        int32_t new_j = st_offset + 1;
                        max_j = std::max(max_j, new_j);
                    } 
                    // ケース2: ノード境界での遷移 (追加)
                    else {
                        int32_t current_k = WavefrontArray::calc_k_from_vk(st_vk);
                        // 次の k = i - next_j = (current_k + offset) - (-1) = current_k + offset + 1
                        int32_t next_k = current_k + st_offset;

                        #pragma GCC unroll 5
                        for (uint8_t code = 1; code <= PatriciaTree::CODE_MAX; ++code) {
                            uint32_t child = _patricia_tree.transition(node_id, code);
                            if (child != 0 && active_counts[child] > 0) {
                                pending_d_buffer[code - 1].push_back_state(child, next_k, 0);
                                has_pending_d = true;
                            }
                        }
                    }
                    idx_d++;
                }

                // --- Deletion (M -> D) ---
                if (k_m_target == min_k) {
                    const int32_t st_offset = wf_array_m.get_offset(idx_m);
                    if (st_offset + 1 < label_len) {
                        int32_t new_j = st_offset + 1;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_m++;
                }

                // 3. 結果の登録
                if (max_j >= -1) {
                    next_wf_array_d.push_back_state(node_id, min_k, max_j);
                }
            }
        }

        idx_m = start_idx_m;  
        size_t idx_i = start_idx_i;  
        
        while (idx_m < end_idx_m || idx_i < end_idx_i) {
            int32_t k_m_target = (idx_m < end_idx_m) ? WavefrontArray::calc_k_from_vk(wf_array_m.get_vk(idx_m)) + 1 : INT32_MAX;
            int32_t k_i_target = (idx_i < end_idx_i) ? WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) + 1 : INT32_MAX;

            int32_t min_k = std::min({k_m_target, k_i_target});
            int32_t max_j = INT32_MIN; 

            if (min_k != INT32_MAX) {
                // --- Insertion (M -> I) ---
                if (k_m_target == min_k) {
                    const uint64_t st_vk = wf_array_m.get_vk(idx_m);
                    const int32_t st_offset = wf_array_m.get_offset(idx_m);
                    int32_t i_pos = WavefrontArray::calc_k_from_vk(st_vk) + st_offset;
                    if (i_pos + 1 < query_length) {
                        int32_t new_j = st_offset;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_m++;
                }

                // --- Insertion (I -> I) ---
                if (k_i_target == min_k) {
                    const uint64_t st_vk = wf_array_i.get_vk(idx_i);
                    const int32_t st_offset = wf_array_i.get_offset(idx_i);
                    int32_t i_pos = WavefrontArray::calc_k_from_vk(st_vk) + st_offset;
                    if (i_pos + 1 < query_length) {
                        int32_t new_j = st_offset;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_i++;
                }

                // 3. 結果の登録
                if (max_j >= -1) {
                    next_wf_array_i.push_back_state(node_id, min_k, max_j);
                }
            }
        }

        if (has_pending_d) {
            #pragma GCC unroll 5
            for (int i = 0; i < PatriciaTree::CODE_MAX; ++i) {
                for (size_t j = 0; j < pending_d_buffer[i].active_size(); ++j) {
                    pending_d.push_back_state(
                        pending_d_buffer[i].get_vk(j),
                        pending_d_buffer[i].get_offset(j)
                    );
                }
                pending_d_buffer[i].clear_logical_size();
            }
            has_pending_d = false;
        }

        start_idx_d = end_idx_d;
        start_idx_m = end_idx_m;
        start_idx_i = end_idx_i;
    }

    // === 保留していたノード境界Deletionのマージ処理 ===
    if (!pending_d.empty()) {
        size_t idx1 = 0;
        size_t idx2 = 0;
        size_t len1 = next_wf_array_d.active_size();
        size_t len2 = pending_d.active_size();

        while (idx1 < len1 || idx2 < len2) {
            if (idx1 == len1) {
                const uint64_t st_vk = pending_d.get_vk(idx2);
                const int32_t st_offset = pending_d.get_offset(idx2);
                merged_wf_array_d.push_back_state(st_vk, st_offset);
                ++idx2;
            } else if (idx2 == len2) {
                // const auto& st = next_wf_array_d[idx1];
                const uint64_t st_vk = next_wf_array_d.get_vk(idx1);
                const int32_t st_offset = next_wf_array_d.get_offset(idx1);
                merged_wf_array_d.push_back_state(st_vk, st_offset);
                ++idx1;
            } else {
                const uint64_t st1_vk = next_wf_array_d.get_vk(idx1);
                const int32_t st1_offset = next_wf_array_d.get_offset(idx1);
                const uint64_t st2_vk = pending_d.get_vk(idx2);
                const int32_t st2_offset = pending_d.get_offset(idx2);

                if (st1_vk < st2_vk) {
                    merged_wf_array_d.push_back_state(st1_vk, st1_offset);
                    ++idx1;
                } else if (st2_vk < st1_vk) {
                    merged_wf_array_d.push_back_state(st2_vk, st2_offset);
                    ++idx2;
                } else {
                    // vk が衝突した場合、offset が大きい（より遠くまで進んでいる）波面を採用
                    const int32_t winner_offset = (st1_offset >= st2_offset) ? st1_offset : st2_offset;
                    merged_wf_array_d.push_back_state(st1_vk, winner_offset);
                    ++idx1;
                    ++idx2;
                }
            }
        }
        // マージ結果を next_wf_array_d に O(1) でスワップ
        next_wf_array_d.swap(merged_wf_array_d);
        merged_wf_array_d.clear_logical_size();
        pending_d.clear_logical_size();
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

    while (start_idx_d < wf_array_dm.active_size() || start_idx_m < wf_array_mm.active_size() || start_idx_i < wf_array_im.active_size()) {
        // 次のnode_idを決定
        uint32_t node_id = UINT32_MAX;
        if (start_idx_d < wf_array_dm.active_size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_dm.get_vk(start_idx_d));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_m < wf_array_mm.active_size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_mm.get_vk(start_idx_m));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_i < wf_array_im.active_size()) {
            uint32_t nid = WavefrontArray::calc_node_id_from_vk(wf_array_im.get_vk(start_idx_i));
            node_id = std::min(node_id, nid);
        }

        size_t end_idx_d = start_idx_d;
        while (end_idx_d < wf_array_dm.active_size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_dm.get_vk(end_idx_d)) == node_id) {
            ++end_idx_d;
        }
        size_t end_idx_m = start_idx_m;
        while (end_idx_m < wf_array_mm.active_size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_mm.get_vk(end_idx_m)) == node_id) {
            ++end_idx_m;
        }
        size_t end_idx_i = start_idx_i;
        while (end_idx_i < wf_array_im.active_size() &&
                WavefrontArray::calc_node_id_from_vk(wf_array_im.get_vk(end_idx_i)) == node_id) {
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
            int32_t k_d_target = (idx_d < end_idx_d) ? WavefrontArray::calc_k_from_vk(wf_array_dm.get_vk(idx_d))     : INT32_MAX;
            int32_t k_m_target = (idx_m < end_idx_m) ? WavefrontArray::calc_k_from_vk(wf_array_mm.get_vk(idx_m))     : INT32_MAX;
            int32_t k_i_target = (idx_i < end_idx_i) ? WavefrontArray::calc_k_from_vk(wf_array_im.get_vk(idx_i))     : INT32_MAX;

            int32_t min_k = std::min({k_d_target, k_m_target, k_i_target});
            int32_t max_j = INT32_MIN; 

            if (min_k != INT32_MAX) {
                if (k_d_target == min_k) {
                    const int32_t st_offset = wf_array_dm.get_offset(idx_d);
                    int32_t new_j = st_offset;
                    max_j = std::max(max_j, new_j);
                    idx_d++;
                }

                if (k_m_target == min_k) {
                    const uint64_t st_vk = wf_array_mm.get_vk(idx_m);
                    const int32_t st_offset = wf_array_mm.get_offset(idx_m);
                    int32_t i_pos = WavefrontArray::calc_k_from_vk(st_vk) + st_offset;
                    if (i_pos + 1 < query_length && st_offset + 1 < label_len) {
                        int32_t new_j = st_offset + 1;
                        max_j = std::max(max_j, new_j);
                    }
                    idx_m++;
                }

                if (k_i_target == min_k) {
                    const int32_t st_offset = wf_array_im.get_offset(idx_i);
                    int32_t new_j = st_offset;
                    max_j = std::max(max_j, new_j);
                    idx_i++;
                }

                if (max_j >= -1) {
                    // if (visited_map_m.update_and_check(WavefrontArray::calc_vk(node_id, min_k), max_j)) {
                        next_wf_array_m.push_back_state(node_id, min_k, max_j);
                    // }
                }
            }
        }
        start_idx_d = end_idx_d;
        start_idx_m = end_idx_m;
        start_idx_i = end_idx_i;
    }
    wf_history_m[(curr_idx + 1) % history_size].swap(next_wf_array_m);
}
