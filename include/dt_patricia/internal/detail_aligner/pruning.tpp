#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// =========================================================
// prune_by_upper_bound
// =========================================================
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::prune_by_upper_bound(
    internal::WavefrontArray &wf_array, const std::vector<uint32_t> &subtree_max_lengths,
    const std::vector<uint32_t> &subtree_min_lengths, int32_t query_length,
    int upper_bound_remain) const
    requires(CostType::is_linear)
{
    size_t write_idx = 0;
    const int32_t max_diff = upper_bound_remain / static_cast<int32_t>(_cost.gap);
    for (size_t i = 0; i < wf_array.active_size(); i++) {
        const uint64_t vk = wf_array.get_vk(i);
        uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vk(vk);
        int32_t k = internal::WavefrontArray::calc_k_from_vk(vk);
        int32_t j_pos = wf_array.get_offset(i);
        int32_t i_pos = k + j_pos;
        int32_t max_remain = subtree_max_lengths[node_id] - (j_pos + 1);
        int32_t min_remain = subtree_min_lengths[node_id] - (j_pos + 1);
        int32_t query_remain = query_length - (i_pos + 1);
        if ((query_remain - max_remain) <= max_diff && (min_remain - query_remain) <= max_diff) {
            if (write_idx != i) {
                wf_array.update_state(write_idx, vk, j_pos);
            }
            write_idx++;
        }
    }
    wf_array.set_size(write_idx);
}

template <AlphabetPolicy Alphabet, typename CostType>
template <bool only_m>
void DTPatricia<Alphabet, CostType>::prune_by_upper_bound(
    internal::WavefrontArray &wf_array_d, internal::WavefrontArray &wf_array_m,
    internal::WavefrontArray &wf_array_i, const std::vector<uint32_t> &subtree_max_lengths,
    const std::vector<uint32_t> &subtree_min_lengths, int32_t query_length,
    int upper_bound_remain) const
    requires(!CostType::is_linear)
{
    const int32_t gap_o = static_cast<int32_t>(_cost.gap_open);
    const int32_t gap_e = static_cast<int32_t>(_cost.gap_extend);

    // M, I, D それぞれに対して判定を行う
    auto prune_logic = [&](internal::WavefrontArray &wf, int state_type) {
        size_t write_idx = 0;
        for (size_t i = 0; i < wf.active_size(); ++i) {
            const uint64_t vk = wf.get_vk(i);
            uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vk(vk);
            int32_t v_k = internal::WavefrontArray::calc_k_from_vk(vk);
            int32_t j_pos = wf.get_offset(i);
            int32_t i_pos = v_k + j_pos;

            int32_t max_rem_t = static_cast<int32_t>(subtree_max_lengths[node_id]) - j_pos;
            int32_t min_rem_t = static_cast<int32_t>(subtree_min_lengths[node_id]) - j_pos;
            int32_t rem_q = query_length - i_pos;

            int32_t lb = 0;
            // クエリがターゲットより長い (Insertionが必要)
            if (rem_q > max_rem_t) {
                int32_t diff = rem_q - max_rem_t;
                // 現在が I 以外なら新たに gap_open が必要
                lb = (state_type == 2 /*STATE_I*/) ? (diff * gap_e) : (gap_o + diff * gap_e);
            }
            // ターゲットがクエリより長い (Deletionが必要)
            else if (min_rem_t > rem_q) {
                int32_t diff = min_rem_t - rem_q;
                // 現在が D 以外なら新たに gap_open が必要
                lb = (state_type == 1 /*STATE_D*/) ? (diff * gap_e) : (gap_o + diff * gap_e);
            }

            if (lb <= upper_bound_remain) {
                if (write_idx != i) {
                    wf.update_state(write_idx, vk, j_pos);
                }
                write_idx++;
            }
        }
        wf.set_size(write_idx);
    };

    prune_logic(wf_array_m, 0);  // Match
    if constexpr (!only_m) {
        prune_logic(wf_array_d, 1);  // Deletion
        prune_logic(wf_array_i, 2);  // Insertion
    }
}

}  // namespace dt_patricia
