#include "PatriciaWFA.hpp"

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
    
    //[NOTE] ここにおいてbufferは常に空かも、検討後に削除するかもしれない
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