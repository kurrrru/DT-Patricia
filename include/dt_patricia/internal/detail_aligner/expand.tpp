#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// ==============================================================================
// Algorithm 3: DT-Patricia Expand - Edit operations (I/D/S)
// ==============================================================================
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::expand(
    const std::string_view query,
    std::vector<internal::WavefrontArray> &wf_history,
    internal::WavefrontArray &next_wf_array,
    int32_t curr_idx,
    size_t history_size,
    const std::vector<uint32_t> &active_counts,
    std::vector<int32_t> &expand_scratch) const requires (CostType::is_linear) {
    next_wf_array.clear_logical_size();
    const int32_t query_length = static_cast<int32_t>(query.length());

    // Sentinel: NULL_OFF < any valid offset; NULL_OFF + 1 is still < 0, so max() naturally rejects it
    static constexpr int32_t NULL_OFF = -1'000'000'000;

    if constexpr (CostType::is_unit) {
        internal::WavefrontArray &wf_array = wf_history[curr_idx];
        next_wf_array.reserve_capacity(wf_array.active_size() * 3 + 10);

        for (size_t idx = 0; idx < wf_array.active_size();) {
            const uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vk(wf_array.get_vk(idx));

            if (active_counts[node_id] == 0) {
                ++idx;
                while (idx < wf_array.active_size() &&
                       internal::WavefrontArray::calc_node_id_from_vk(wf_array.get_vk(idx)) == node_id)
                    ++idx;
                continue;
            }

            const size_t start_idx = idx;
            size_t end_idx = idx + 1;
            while (end_idx < wf_array.active_size() &&
                   internal::WavefrontArray::calc_node_id_from_vk(wf_array.get_vk(end_idx)) == node_id)
                ++end_idx;

            const int32_t label_len = static_cast<int32_t>(_patricia_tree.get_label_length(node_id));
            const int32_t k_lo = internal::WavefrontArray::calc_k_from_vk(wf_array.get_vk(start_idx));
            const int32_t k_hi = internal::WavefrontArray::calc_k_from_vk(wf_array.get_vk(end_idx - 1));

            // scratch indexed by (k - k_base), k_base = k_lo - 2
            // [0]          = k_lo-2  sentinel (always NULL)
            // [1]          = k_lo-1  left pad (NULL after scatter)
            // [2..range+1] = k_lo..k_hi
            // [range+2]    = k_hi+1  right pad (NULL after scatter)
            // [range+3]    = k_hi+2  sentinel for del read at k=k_hi+1
            const int32_t k_base = k_lo - 2;
            const int32_t scratch_sz = k_hi - k_lo + 5;

            if ((int32_t)expand_scratch.size() < scratch_sz)
                expand_scratch.assign(scratch_sz, NULL_OFF);
            else
                std::fill(expand_scratch.begin(), expand_scratch.begin() + scratch_sz, NULL_OFF);

            for (size_t i = start_idx; i < end_idx; ++i) {
                int32_t k = internal::WavefrontArray::calc_k_from_vk(wf_array.get_vk(i));
                int32_t j = wf_array.get_offset(i);
                int32_t& slot = expand_scratch[k - k_base];
                if (j > slot) slot = j;  // take max on duplicate k (shouldn't happen but safe)
            }

            // Dense loop over output range [k_lo-1, k_hi+1]
            for (int32_t k = k_lo - 1; k <= k_hi + 1; ++k) {
                const int32_t s = k - k_base;  // s in [1, k_hi-k_lo+3]

                const int32_t src_d = expand_scratch[s + 1];
                const int32_t del_j = (src_d > NULL_OFF && src_d + 1 < label_len)
                                      ? src_d + 1 : INT32_MIN;

                const int32_t src_s = expand_scratch[s];
                const int32_t sub_j = (src_s > NULL_OFF && src_s + 1 < label_len
                                       && k + src_s + 1 < query_length)
                                      ? src_s + 1 : INT32_MIN;

                const int32_t src_i = expand_scratch[s - 1];
                const int32_t ins_j = (src_i > NULL_OFF && k - 1 + src_i + 1 < query_length)
                                      ? src_i : INT32_MIN;

                const int32_t max_j = std::max(del_j, std::max(sub_j, ins_j));
                if (max_j >= -1)
                    next_wf_array.push_back_unchecked(internal::WavefrontArray::calc_vk(node_id, k), max_j);
            }

            idx = end_idx;
        }
        wf_history[internal::add_mod(curr_idx, 1, history_size)].swap(next_wf_array);

    } else {
        // Linear gap cost: 3 separate source arrays — 3-pointer merge (fill overhead of scatter
        // dominates for large k when density is low across 3 independent arrays)
        size_t wf_history_idx_d = internal::sub_mod(curr_idx + 1, _cost.gap, history_size);
        size_t wf_history_idx_s = internal::sub_mod(curr_idx + 1, _cost.mismatch, history_size);
        size_t wf_history_idx_i = internal::sub_mod(curr_idx + 1, _cost.gap, history_size);

        internal::WavefrontArray &wf_array_d = wf_history[wf_history_idx_d];
        internal::WavefrontArray &wf_array_s = wf_history[wf_history_idx_s];
        internal::WavefrontArray &wf_array_i = wf_history[wf_history_idx_i];

        next_wf_array.reserve_capacity(
            (wf_array_d.active_size() + wf_array_s.active_size() + wf_array_i.active_size()) * 3 + 10);

        size_t start_idx_d = 0, start_idx_s = 0, start_idx_i = 0;

        while (start_idx_d < wf_array_d.active_size() ||
               start_idx_s < wf_array_s.active_size() ||
               start_idx_i < wf_array_i.active_size()) {

            uint32_t node_id = UINT32_MAX;
            if (start_idx_d < wf_array_d.active_size())
                node_id = std::min(node_id, internal::WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(start_idx_d)));
            if (start_idx_s < wf_array_s.active_size())
                node_id = std::min(node_id, internal::WavefrontArray::calc_node_id_from_vk(wf_array_s.get_vk(start_idx_s)));
            if (start_idx_i < wf_array_i.active_size())
                node_id = std::min(node_id, internal::WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(start_idx_i)));

            size_t end_idx_d = start_idx_d;
            while (end_idx_d < wf_array_d.active_size() &&
                   internal::WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(end_idx_d)) == node_id)
                ++end_idx_d;
            size_t end_idx_s = start_idx_s;
            while (end_idx_s < wf_array_s.active_size() &&
                   internal::WavefrontArray::calc_node_id_from_vk(wf_array_s.get_vk(end_idx_s)) == node_id)
                ++end_idx_s;
            size_t end_idx_i = start_idx_i;
            while (end_idx_i < wf_array_i.active_size() &&
                   internal::WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(end_idx_i)) == node_id)
                ++end_idx_i;

            if (active_counts[node_id] == 0) {
                start_idx_d = end_idx_d;
                start_idx_s = end_idx_s;
                start_idx_i = end_idx_i;
                continue;
            }

            const int32_t label_len = static_cast<int32_t>(_patricia_tree.get_label_length(node_id));

            size_t idx_d = start_idx_d;
            size_t idx_s = start_idx_s;
            size_t idx_i = start_idx_i;

            int32_t k_d_raw = (idx_d < end_idx_d) ? internal::WavefrontArray::calc_k_from_vk(wf_array_d.get_vk(idx_d)) : 0;
            int32_t k_s_raw = (idx_s < end_idx_s) ? internal::WavefrontArray::calc_k_from_vk(wf_array_s.get_vk(idx_s)) : 0;
            int32_t k_i_raw = (idx_i < end_idx_i) ? internal::WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) : 0;

            while (idx_d < end_idx_d || idx_s < end_idx_s || idx_i < end_idx_i) {
                const int32_t k_d_target = (idx_d < end_idx_d) ? k_d_raw - 1 : INT32_MAX;
                const int32_t k_s_target = (idx_s < end_idx_s) ? k_s_raw     : INT32_MAX;
                const int32_t k_i_target = (idx_i < end_idx_i) ? k_i_raw + 1 : INT32_MAX;

                const int32_t min_k = std::min(k_d_target, std::min(k_s_target, k_i_target));
                int32_t max_j = INT32_MIN;

                if (k_d_target == min_k) {
                    const int32_t st_offset = wf_array_d.get_offset(idx_d);
                    if (st_offset + 1 < label_len) max_j = st_offset + 1;
                    ++idx_d;
                    k_d_raw = (idx_d < end_idx_d) ? internal::WavefrontArray::calc_k_from_vk(wf_array_d.get_vk(idx_d)) : 0;
                }

                if (k_s_target == min_k) {
                    const int32_t st_offset = wf_array_s.get_offset(idx_s);
                    const int32_t i_pos = k_s_raw + st_offset;
                    if (i_pos + 1 < query_length && st_offset + 1 < label_len) {
                        const int32_t new_j = st_offset + 1;
                        if (new_j > max_j) max_j = new_j;
                    }
                    ++idx_s;
                    k_s_raw = (idx_s < end_idx_s) ? internal::WavefrontArray::calc_k_from_vk(wf_array_s.get_vk(idx_s)) : 0;
                }

                if (k_i_target == min_k) {
                    const int32_t st_offset = wf_array_i.get_offset(idx_i);
                    const int32_t i_pos = k_i_raw + st_offset;
                    if (i_pos + 1 < query_length && st_offset > max_j) max_j = st_offset;
                    ++idx_i;
                    k_i_raw = (idx_i < end_idx_i) ? internal::WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) : 0;
                }

                if (max_j >= -1)
                    next_wf_array.push_back_unchecked(internal::WavefrontArray::calc_vk(node_id, min_k), max_j);
            }

            start_idx_d = end_idx_d;
            start_idx_s = end_idx_s;
            start_idx_i = end_idx_i;
        }
        wf_history[internal::add_mod(curr_idx, 1, history_size)].swap(next_wf_array);
    }
}

// ==============================================================================
// Algorithm 3: DT-Patricia Expand - Edit operations (I/D/S) for affine gap cost
// ==============================================================================
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::expand(
    const std::string_view query,
    std::vector<internal::WavefrontArray> &wf_history_d,
    std::vector<internal::WavefrontArray> &wf_history_m,
    std::vector<internal::WavefrontArray> &wf_history_i,
    internal::WavefrontArray &next_wf_array_d,
    internal::WavefrontArray &next_wf_array_m,
    internal::WavefrontArray &next_wf_array_i,
    int32_t curr_idx,
    size_t history_size,
    const std::vector<uint32_t> &active_counts,
    std::array<internal::WavefrontArray, PatriciaTree<Alphabet>::CODE_MAX> &pending_d_buffer,
    internal::WavefrontArray &pending_d,
    internal::WavefrontArray &merged_wf_array_d,
    std::vector<int32_t> &expand_scratch) const requires (!CostType::is_linear) {
    (void)expand_scratch;
    next_wf_array_d.clear_logical_size();
    next_wf_array_m.clear_logical_size();
    next_wf_array_i.clear_logical_size();
    const int32_t query_length = static_cast<int32_t>(query.length());

    const uint32_t next_idx = internal::add_mod(curr_idx, 1, history_size);
    // current_score + 1 に対応する internal::WavefrontArray を構築する
    size_t wf_history_idx_d = internal::sub_mod(next_idx, _cost.gap_extend, history_size);
    size_t wf_history_idx_m = internal::sub_mod(next_idx, _cost.gap_open + _cost.gap_extend, history_size);
    size_t wf_history_idx_i = internal::sub_mod(next_idx, _cost.gap_extend, history_size);

    internal::WavefrontArray &wf_array_d = wf_history_d[wf_history_idx_d];
    internal::WavefrontArray &wf_array_m = wf_history_m[wf_history_idx_m];
    internal::WavefrontArray &wf_array_i = wf_history_i[wf_history_idx_i];

    // Pre-reserve capacity for deletion and insertion outputs
    next_wf_array_d.reserve_capacity(
        (wf_array_d.active_size() + wf_array_m.active_size()) * 3 + 10);
    next_wf_array_i.reserve_capacity(
        (wf_array_m.active_size() + wf_array_i.active_size()) * 3 + 10);

    size_t start_idx_d = 0;
    size_t start_idx_m = 0;
    size_t start_idx_i = 0;

    bool has_pending_d = false;

    while (start_idx_d < wf_array_d.active_size() || start_idx_m < wf_array_m.active_size() || start_idx_i < wf_array_i.active_size()) {
        // 次のnode_idを決定
        uint32_t node_id = UINT32_MAX;
        if (start_idx_d < wf_array_d.active_size()) {
            uint32_t nid = internal::WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(start_idx_d));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_m < wf_array_m.active_size()) {
            uint32_t nid = internal::WavefrontArray::calc_node_id_from_vk(wf_array_m.get_vk(start_idx_m));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_i < wf_array_i.active_size()) {
            uint32_t nid = internal::WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(start_idx_i));
            node_id = std::min(node_id, nid);
        }

        // node_idが等しい区間を見つける
        size_t end_idx_d = start_idx_d;
        while (end_idx_d < wf_array_d.active_size() &&
                internal::WavefrontArray::calc_node_id_from_vk(wf_array_d.get_vk(end_idx_d)) == node_id) {
            ++end_idx_d;
        }
        size_t end_idx_m = start_idx_m;
        while (end_idx_m < wf_array_m.active_size() &&
                internal::WavefrontArray::calc_node_id_from_vk(wf_array_m.get_vk(end_idx_m)) == node_id) {
            ++end_idx_m;
        }
        size_t end_idx_i = start_idx_i;
        while (end_idx_i < wf_array_i.active_size() &&
                internal::WavefrontArray::calc_node_id_from_vk(wf_array_i.get_vk(end_idx_i)) == node_id) {
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

        // Cache k values for D-stream merge
        int32_t k_d_raw = (idx_d < end_idx_d) ? internal::WavefrontArray::calc_k_from_vk(wf_array_d.get_vk(idx_d)) : 0;
        int32_t k_m_d_raw = (idx_m < end_idx_m) ? internal::WavefrontArray::calc_k_from_vk(wf_array_m.get_vk(idx_m)) : 0;

        while (idx_d < end_idx_d || idx_m < end_idx_m) {
            const int32_t k_d_target = (idx_d < end_idx_d) ? k_d_raw - 1 : INT32_MAX;
            const int32_t k_m_target = (idx_m < end_idx_m) ? k_m_d_raw - 1 : INT32_MAX;

            const int32_t min_k = std::min(k_d_target, k_m_target);
            int32_t max_j = INT32_MIN;

            // --- Deletion (D -> D) ---
            if (k_d_target == min_k) {
                const uint64_t st_vk = wf_array_d.get_vk(idx_d);
                const int32_t st_offset = wf_array_d.get_offset(idx_d);
                // ケース1: ノード内での伸長
                if (st_offset + 1 < label_len) {
                    const int32_t new_j = st_offset + 1;
                    max_j = new_j;
                }
                // ケース2: ノード境界での遷移 (追加)
                else {
                    const int32_t current_k = internal::WavefrontArray::calc_k_from_vk(st_vk);
                    const int32_t next_k = current_k + st_offset;

                    for (uint8_t code = 1; code <= tree_type::CODE_MAX; ++code) {
                        uint32_t child = _patricia_tree.transition(node_id, code);
                        if (child != 0 && active_counts[child] > 0) {
                            pending_d_buffer[code - 1].push_back_state(child, next_k, 0);
                            has_pending_d = true;
                        }
                    }
                }
                ++idx_d;
                k_d_raw = (idx_d < end_idx_d) ? internal::WavefrontArray::calc_k_from_vk(wf_array_d.get_vk(idx_d)) : 0;
            }

            // --- Deletion (M -> D) ---
            if (k_m_target == min_k) {
                const int32_t st_offset = wf_array_m.get_offset(idx_m);
                if (st_offset + 1 < label_len) {
                    const int32_t new_j = st_offset + 1;
                    if (new_j > max_j) max_j = new_j;
                }
                ++idx_m;
                k_m_d_raw = (idx_m < end_idx_m) ? internal::WavefrontArray::calc_k_from_vk(wf_array_m.get_vk(idx_m)) : 0;
            }

            // 3. 結果の登録
            if (max_j >= -1) {
                next_wf_array_d.push_back_unchecked(internal::WavefrontArray::calc_vk(node_id, min_k), max_j);
            }
        }

        idx_m = start_idx_m;
        size_t idx_i = start_idx_i;

        // Cache k values for I-stream merge
        int32_t k_m_i_raw = (idx_m < end_idx_m) ? internal::WavefrontArray::calc_k_from_vk(wf_array_m.get_vk(idx_m)) : 0;
        int32_t k_i_raw = (idx_i < end_idx_i) ? internal::WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) : 0;

        while (idx_m < end_idx_m || idx_i < end_idx_i) {
            const int32_t k_m_target = (idx_m < end_idx_m) ? k_m_i_raw + 1 : INT32_MAX;
            const int32_t k_i_target = (idx_i < end_idx_i) ? k_i_raw + 1 : INT32_MAX;

            const int32_t min_k = std::min(k_m_target, k_i_target);
            int32_t max_j = INT32_MIN;

            // --- Insertion (M -> I) ---
            if (k_m_target == min_k) {
                const int32_t st_offset = wf_array_m.get_offset(idx_m);
                const int32_t i_pos = k_m_i_raw + st_offset;
                if (i_pos + 1 < query_length) {
                    const int32_t new_j = st_offset;
                    max_j = new_j;
                }
                ++idx_m;
                k_m_i_raw = (idx_m < end_idx_m) ? internal::WavefrontArray::calc_k_from_vk(wf_array_m.get_vk(idx_m)) : 0;
            }

            // --- Insertion (I -> I) ---
            if (k_i_target == min_k) {
                const int32_t st_offset = wf_array_i.get_offset(idx_i);
                const int32_t i_pos = k_i_raw + st_offset;
                if (i_pos + 1 < query_length) {
                    const int32_t new_j = st_offset;
                    if (new_j > max_j) max_j = new_j;
                }
                ++idx_i;
                k_i_raw = (idx_i < end_idx_i) ? internal::WavefrontArray::calc_k_from_vk(wf_array_i.get_vk(idx_i)) : 0;
            }

            // 3. 結果の登録
            if (max_j >= -1) {
                next_wf_array_i.push_back_unchecked(internal::WavefrontArray::calc_vk(node_id, min_k), max_j);
            }
        }

        if (has_pending_d) {
            #pragma GCC unroll 5
            for (int i = 0; i < tree_type::CODE_MAX; ++i) {
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


    wf_history_d[next_idx].swap(next_wf_array_d);
    wf_history_i[next_idx].swap(next_wf_array_i);

    // match/mismatch のみを見る
    wf_history_idx_d = next_idx;
    wf_history_idx_m = internal::sub_mod(next_idx, _cost.mismatch, history_size);
    wf_history_idx_i = next_idx;

    internal::WavefrontArray &wf_array_dm = wf_history_d[wf_history_idx_d];
    internal::WavefrontArray &wf_array_mm = wf_history_m[wf_history_idx_m];
    internal::WavefrontArray &wf_array_im = wf_history_i[wf_history_idx_i];

    // Pre-reserve capacity for the M-update pass
    next_wf_array_m.reserve_capacity(
        (wf_array_dm.active_size() + wf_array_mm.active_size() + wf_array_im.active_size()) * 3 + 10);

    start_idx_d = 0;
    start_idx_m = 0;
    start_idx_i = 0;

    while (start_idx_d < wf_array_dm.active_size() || start_idx_m < wf_array_mm.active_size() || start_idx_i < wf_array_im.active_size()) {
        // 次のnode_idを決定
        uint32_t node_id = UINT32_MAX;
        if (start_idx_d < wf_array_dm.active_size()) {
            uint32_t nid = internal::WavefrontArray::calc_node_id_from_vk(wf_array_dm.get_vk(start_idx_d));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_m < wf_array_mm.active_size()) {
            uint32_t nid = internal::WavefrontArray::calc_node_id_from_vk(wf_array_mm.get_vk(start_idx_m));
            node_id = std::min(node_id, nid);
        }
        if (start_idx_i < wf_array_im.active_size()) {
            uint32_t nid = internal::WavefrontArray::calc_node_id_from_vk(wf_array_im.get_vk(start_idx_i));
            node_id = std::min(node_id, nid);
        }

        size_t end_idx_d = start_idx_d;
        while (end_idx_d < wf_array_dm.active_size() &&
                internal::WavefrontArray::calc_node_id_from_vk(wf_array_dm.get_vk(end_idx_d)) == node_id) {
            ++end_idx_d;
        }
        size_t end_idx_m = start_idx_m;
        while (end_idx_m < wf_array_mm.active_size() &&
                internal::WavefrontArray::calc_node_id_from_vk(wf_array_mm.get_vk(end_idx_m)) == node_id) {
            ++end_idx_m;
        }
        size_t end_idx_i = start_idx_i;
        while (end_idx_i < wf_array_im.active_size() &&
                internal::WavefrontArray::calc_node_id_from_vk(wf_array_im.get_vk(end_idx_i)) == node_id) {
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

        // Cache k values for M-update pass
        int32_t k_d_raw = (idx_d < end_idx_d) ? internal::WavefrontArray::calc_k_from_vk(wf_array_dm.get_vk(idx_d)) : 0;
        int32_t k_m_raw = (idx_m < end_idx_m) ? internal::WavefrontArray::calc_k_from_vk(wf_array_mm.get_vk(idx_m)) : 0;
        int32_t k_i_raw = (idx_i < end_idx_i) ? internal::WavefrontArray::calc_k_from_vk(wf_array_im.get_vk(idx_i)) : 0;

        while (idx_d < end_idx_d || idx_m < end_idx_m || idx_i < end_idx_i) {
            const int32_t k_d_target = (idx_d < end_idx_d) ? k_d_raw : INT32_MAX;
            const int32_t k_m_target = (idx_m < end_idx_m) ? k_m_raw : INT32_MAX;
            const int32_t k_i_target = (idx_i < end_idx_i) ? k_i_raw : INT32_MAX;

            const int32_t min_k = std::min(k_d_target, std::min(k_m_target, k_i_target));
            int32_t max_j = INT32_MIN;

            if (k_d_target == min_k) {
                const int32_t st_offset = wf_array_dm.get_offset(idx_d);
                const int32_t new_j = st_offset;
                max_j = new_j;
                ++idx_d;
                k_d_raw = (idx_d < end_idx_d) ? internal::WavefrontArray::calc_k_from_vk(wf_array_dm.get_vk(idx_d)) : 0;
            }

            if (k_m_target == min_k) {
                const int32_t st_offset = wf_array_mm.get_offset(idx_m);
                const int32_t i_pos = k_m_raw + st_offset;
                if (i_pos + 1 < query_length && st_offset + 1 < label_len) {
                    const int32_t new_j = st_offset + 1;
                    if (new_j > max_j) max_j = new_j;
                }
                ++idx_m;
                k_m_raw = (idx_m < end_idx_m) ? internal::WavefrontArray::calc_k_from_vk(wf_array_mm.get_vk(idx_m)) : 0;
            }

            if (k_i_target == min_k) {
                const int32_t st_offset = wf_array_im.get_offset(idx_i);
                const int32_t new_j = st_offset;
                if (new_j > max_j) max_j = new_j;
                ++idx_i;
                k_i_raw = (idx_i < end_idx_i) ? internal::WavefrontArray::calc_k_from_vk(wf_array_im.get_vk(idx_i)) : 0;
            }

            if (max_j >= -1) {
                next_wf_array_m.push_back_unchecked(internal::WavefrontArray::calc_vk(node_id, min_k), max_j);
            }
        }
        start_idx_d = end_idx_d;
        start_idx_m = end_idx_m;
        start_idx_i = end_idx_i;
    }
    wf_history_m[next_idx].swap(next_wf_array_m);
}

}  // namespace dt_patricia
