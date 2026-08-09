#include <cassert>

#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// =========================================================
// Algorithm 2: DT-Patricia Extend - Exact match extension
// =========================================================
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::extend(
    const std::string_view query, internal::WavefrontArray &wf_array,
    internal::WavefrontArray &next_wf_array, internal::WavefrontArray &child_wf_array,
    std::array<internal::WavefrontArray, PatriciaTree<Alphabet>::CODE_MAX> &buffer,
    const std::vector<uint32_t> &active_counts, internal::ReachedOffsetTable &reached,
    int32_t current_score) const {
    next_wf_array.clear_logical_size();
    child_wf_array.clear_logical_size();

    bool buffer_used = false;

    const int32_t query_length = static_cast<int32_t>(query.length());

    size_t wf_array_idx = 0;
    size_t child_idx = 0;
    uint32_t last_node_id = UINT32_MAX;

    // Label caching to avoid repeated get_label calls for the same node_id
    std::string_view cached_label;
    int32_t cached_label_len = 0;
    uint32_t cached_label_node_id = UINT32_MAX;

    while (wf_array_idx < wf_array.active_size() || child_idx < child_wf_array.active_size() ||
           buffer_used) {
        bool wf_array_idx_increment = false;
        bool child_idx_increment = false;

        // === ステップ1: 次に処理する状態を選択 ===
        uint64_t current_vk;
        int32_t current_offset;
        if (child_idx >= child_wf_array.active_size() && wf_array_idx >= wf_array.active_size()) {
            // bufferのみ残っている場合
            for (const auto &buf : buffer) {
                for (size_t i = 0; i < buf.active_size(); ++i) {
                    child_wf_array.push_back_state(buf.get_vk(i), buf.get_offset(i));
                }
            }
            for (auto &buf : buffer) {
                buf.clear_logical_size();
            }
            buffer_used = false;
            continue;
        }

        if (child_idx >= child_wf_array.active_size()) {
            current_vk = wf_array.get_vk(wf_array_idx);
            current_offset = wf_array.get_offset(wf_array_idx);
            wf_array_idx_increment = true;
        } else if (wf_array_idx >= wf_array.active_size()) {
            current_vk = child_wf_array.get_vk(child_idx);
            current_offset = child_wf_array.get_offset(child_idx);
            child_idx_increment = true;
        } else {
            const uint64_t wf_vk = wf_array.get_vk(wf_array_idx);
            const uint64_t child_vk = child_wf_array.get_vk(child_idx);

            if (wf_vk < child_vk) {
                current_vk = wf_vk;
                current_offset = wf_array.get_offset(wf_array_idx);
                wf_array_idx_increment = true;
            } else if (wf_vk > child_vk) {
                current_vk = child_vk;
                current_offset = child_wf_array.get_offset(child_idx);
                child_idx_increment = true;
            } else {
                // vk が同じ → offset が大きい方を採用
                current_vk = wf_vk;
                current_offset = std::max(wf_array.get_offset(wf_array_idx),
                                          child_wf_array.get_offset(child_idx));
                wf_array_idx_increment = true;
                child_idx_increment = true;
            }
        }

        // === ステップ2: node_idが変わったらbufferを処理 ===
        uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vk(current_vk);
        if (node_id != last_node_id && buffer_used) {
            // bufferはソート済み
            // (BFS順により child_wf_arrayの末尾 < bufferの最小値 が保証される)
            for (const auto &buf : buffer) {
                for (size_t i = 0; i < buf.active_size(); ++i) {
                    child_wf_array.push_back_state(buf.get_vk(i), buf.get_offset(i));
                }
            }
            for (auto &buf : buffer) {
                buf.clear_logical_size();
            }
            buffer_used = false;

            // これによってchild_wf_arrayに追加された値の方がwf_arrayの値より小さいかもしれない
            // child_idx_incrementが立っていない = child_idxの値はチェックされていない可能性がある
            if (!child_idx_increment) {
                const uint64_t child_vk = child_wf_array.get_vk(child_idx);
                const int32_t child_offset = child_wf_array.get_offset(child_idx);
                if (child_vk < current_vk) {
                    current_vk = child_vk;
                    current_offset = child_offset;
                    child_idx_increment = true;
                    wf_array_idx_increment = false;
                } else if (child_vk == current_vk) {
                    if (child_offset > current_offset) {
                        current_vk = child_vk;
                        current_offset = child_offset;
                    }
                    child_idx_increment = true;
                    wf_array_idx_increment = true;
                }
                node_id = internal::WavefrontArray::calc_node_id_from_vk(current_vk);
            }
        }
        last_node_id = node_id;

        // === ステップ3: 支配判定（入口版）===
        int32_t k = internal::WavefrontArray::calc_k_from_vk(current_vk);
        int32_t j = current_offset;

        if constexpr (internal::DOMINANCE_AT_EXTEND_ENTRY) {
            const bool dominated = reached.dominated(node_id, k, j, current_score);
            if (dominated) {
                if (wf_array_idx_increment) {
                    ++wf_array_idx;
                }
                if (child_idx_increment) {
                    ++child_idx;
                }
                continue;
            }
        }

        // === ステップ4: Extension 処理 ===
        int32_t i = k + j;

        // Cache label lookup: only call get_label when node_id changes
        if (node_id != cached_label_node_id) {
            cached_label = _patricia_tree.get_label(node_id);
            cached_label_len = static_cast<int32_t>(cached_label.length());
            cached_label_node_id = node_id;
        }
        const std::string_view &label = cached_label;
        const int32_t label_len = cached_label_len;

        const int32_t max_lcp_len = std::min(query_length - (i + 1), label_len - (j + 1));
        // exact match extension
        int32_t lcp_len = static_cast<int32_t>(
            internal::fast_lcp(query.data() + i + 1, label.data() + j + 1, max_lcp_len));
        i += lcp_len;
        j += lcp_len;

        if constexpr (internal::DOMINANCE_AT_EXTEND_ENTRY) {
            reached.record(node_id, k, j, current_score);  // 伸長後の到達点で更新する
        }

        // === ステップ5: 子ノードへの遷移 or next_wf_array への追加 ===
        if (j + 1 == label_len) {
            bool already_expanded = false;
            if constexpr (!internal::DOMINANCE_AT_EXTEND_ENTRY) {
                already_expanded = reached.dominated(node_id, k, j, current_score);
                if (!already_expanded) {
                    reached.record(node_id, k, j, current_score);
                }
            }

            // ノード終端に到達 → 子ノードをbufferに追加
            for (uint8_t code = 1; !already_expanded && code <= tree_type::CODE_MAX; ++code) {
                uint32_t child = _patricia_tree.transition(node_id, code);
                if (child != 0 && active_counts[child] > 0) {
                    int32_t new_k = (i + 1) - 0;
                    buffer[code - 1].push_back_state(child, new_k, -1);
                    buffer_used = true;
                }
            }

            if (_patricia_tree.is_terminal(node_id)) {
                next_wf_array.push_back_state(node_id, k, j);
            }
        } else {
            next_wf_array.push_back_state(node_id, k, j);
        }
        if (wf_array_idx_increment) {
            ++wf_array_idx;
        }
        if (child_idx_increment) {
            ++child_idx;
        }
    }
    // ここではbufferは空になっている(空になっていないとwhileループが継続するため)
#ifndef NDEBUG
    for (const auto &buf : buffer) {
        assert(buf.empty());
    }
#endif
    wf_array.swap(next_wf_array);
}

}  // namespace dt_patricia
