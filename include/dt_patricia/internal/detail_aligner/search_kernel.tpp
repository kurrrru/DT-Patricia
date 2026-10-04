#include <cassert>

#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// =========================================================
// テンプレート関数の実装
// =========================================================
template <AlphabetPolicy Alphabet, typename CostType>
template <typename StopPredicate>
std::vector<AlignmentResult> DTPatricia<Alphabet, CostType>::search_kernel(
    const std::string &query, StopPredicate stop_predicate, int upper_bound)
    requires(CostType::is_linear)
{
#ifndef NDEBUG
    internal::assert_no_null_bytes(query);
#endif
    std::vector<AlignmentResult> results;

    const std::vector<uint32_t> &subtree_max_lengths = _patricia_tree.get_subtree_max_lengths();
    const std::vector<uint32_t> &subtree_min_lengths = _patricia_tree.get_subtree_min_lengths();

    if (_patricia_tree.empty()) {
        return results;
    }

    // ここから先で書き換えた状態（c, _reported, reached）は、抜けるときに必ず戻す
    SearchStateGuard state_guard(*this, results);

    // =========================================================
    // クエリのパディング処理
    // =========================================================
    std::string padded_query_str = query;
    canonicalize_inplace<Alphabet>(padded_query_str.data(), query.length());
    padded_query_str.append(tree_type::SIMD_PADDING_SIZE, '\0');
    std::string_view padded_query(padded_query_str.data(), query.length());
    const int32_t query_length = static_cast<int32_t>(padded_query.length());

    uint32_t history_size;
    if constexpr (CostType::is_unit) {  // Simple Edit Distances
        history_size = 2;
    } else {  // Linear
        history_size = std::max(_cost.mismatch, _cost.gap) + 1;
    }

    std::vector<internal::WavefrontArray> wf_history(history_size);
    internal::WavefrontArray next_wf_array;
    internal::WavefrontArray child_wf_array;
    std::array<internal::WavefrontArray, tree_type::CODE_MAX> buffer;
    std::vector<int32_t> expand_scratch;
    std::vector<int32_t> expand_maxj;

    // 再処理の抑止に使う表 (_reached)。状態管理のオーバーヘッドがあるので、探索が育って元が取れる
    // 見込みが立つまで有効化しない。
    size_t states_seen = 0;
    const size_t reached_enable_threshold = _patricia_tree.node_count();

    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    const uint32_t root = _patricia_tree.root_id();
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

        if (!_reached.enabled()) {
            states_seen += curr_wf.active_size();
            if (states_seen > reached_enable_threshold) {
                _reached.enable();
            }
        }

        // Algorithm 2: DT-Patricia Extend
        extend(padded_query, curr_wf, next_wf_array, child_wf_array, buffer, _active_counts,
               _reached, current_score);

        if (upper_bound >= 0) {
            prune_by_upper_bound(curr_wf, subtree_max_lengths, subtree_min_lengths, query_length,
                                 upper_bound - current_score);
        }

        // 終端チェック: クエリ全体が処理されたノードを探す。
        // c はこの場で更新する。c を読む expand はこのあとに実行されるので、
        // スコアごとにまとめて反映するのと同じ値を読む。
        for (size_t idx = 0; idx < curr_wf.active_size(); ++idx) {
            uint64_t curr_vd = curr_wf.get_vd(idx);
            uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vd(curr_vd);
            int32_t diag = internal::WavefrontArray::calc_diag_from_vd(curr_vd);
            int32_t j = curr_wf.get_offset(idx);
            int32_t i = diag + j;

            if (i + 1 == query_length &&
                j + 1 == static_cast<int32_t>(_patricia_tree.get_label_length(node_id))) {
                report_strings_at(node_id, current_score, results);
            }
        }

        // 停止条件チェック
        if (stop_predicate(current_score, results)) {
            break;
        }

        // Algorithm 3: DT-Patricia Expand
        expand(padded_query, wf_history, next_wf_array, curr_idx, history_size, _active_counts,
               expand_scratch, expand_maxj);

        uint32_t next_idx = internal::increment_mod(curr_idx, history_size);
        if (upper_bound >= 0) {
            prune_by_upper_bound(wf_history[next_idx], subtree_max_lengths, subtree_min_lengths,
                                 query_length, upper_bound - (current_score + 1));
        }

        ++current_score;
        curr_idx = next_idx;
        next_idx = internal::increment_mod(next_idx, history_size);

        // 以降のループで使うことはないので履歴をリセット
        wf_history[next_idx].clear_logical_size();
    }
    state_guard.restore();
    return results;
}

template <AlphabetPolicy Alphabet, typename CostType>
template <typename StopPredicate>
std::vector<AlignmentResult> DTPatricia<Alphabet, CostType>::search_kernel(
    const std::string &query, StopPredicate stop_predicate, int upper_bound)
    requires(!CostType::is_linear)
{
#ifndef NDEBUG
    internal::assert_no_null_bytes(query);
#endif
    std::vector<AlignmentResult> results;

    if (_patricia_tree.empty()) {
        return results;
    }

    // ここから先で書き換えた状態（c, _reported, reached）は、抜けるときに必ず戻す
    SearchStateGuard state_guard(*this, results);

    // =========================================================
    // クエリのパディング処理
    // =========================================================
    std::string padded_query_str = query;
    canonicalize_inplace<Alphabet>(padded_query_str.data(), query.length());
    padded_query_str.append(tree_type::SIMD_PADDING_SIZE, '\0');
    std::string_view padded_query(padded_query_str.data(), query.length());
    const int32_t query_length = static_cast<int32_t>(padded_query.length());

    uint32_t history_size =
        std::max({_cost.mismatch, _cost.gap_open + _cost.gap_extend, _cost.gap_extend}) + 1;

    std::vector<internal::WavefrontArray> wf_history_d(history_size);
    std::vector<internal::WavefrontArray> wf_history_m(history_size);
    std::vector<internal::WavefrontArray> wf_history_i(history_size);
    internal::WavefrontArray next_wf_array_d;
    internal::WavefrontArray next_wf_array_m;
    internal::WavefrontArray next_wf_array_i;
    internal::WavefrontArray child_wf_array;
    std::array<internal::WavefrontArray, tree_type::CODE_MAX> buffer;

    internal::WavefrontArray pending_d;
    internal::WavefrontArray merged_wf_array_d;
    std::vector<int32_t> expand_scratch;

    // _reached_d: D 層の子生成（expand の pending_d 経由）の重複を潰す表。
    // _reached: M 層の子生成（extend 経由）の重複を潰す表。
    size_t states_seen = 0;
    const size_t reached_enable_threshold = _patricia_tree.node_count();

    // 初期状態: ルートノードから開始 (i=-1, j=-1, diagonal=0)
    const uint32_t root = _patricia_tree.root_id();
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

        // reachedとreached_dは必ず同時に有効化される
        assert(_reached.enabled() == _reached_d.enabled());
        if (!_reached_d.enabled()) {
            states_seen += curr_wf_m.active_size();
            if (states_seen > reached_enable_threshold) {
                _reached_d.enable();
                _reached.enable();
            }
        }

        // Algorithm 2: DT-Patricia Extend
        extend(padded_query, curr_wf_m, next_wf_array_m, child_wf_array, buffer, _active_counts,
               _reached, current_score);
        if (upper_bound >= 0) {
            prune_by_upper_bound<true>(next_wf_array_d, curr_wf_m, next_wf_array_i,
                                       _patricia_tree.get_subtree_max_lengths(),
                                       _patricia_tree.get_subtree_min_lengths(), query_length,
                                       upper_bound - current_score);
        }

        // 終端チェック: クエリ全体が処理されたノードを探す。
        // c はこの場で更新する。c を読む expand はこのあとに実行されるので、
        // スコアごとにまとめて反映するのと同じ値を読む。
        for (size_t idx = 0; idx < curr_wf_m.active_size(); ++idx) {
            uint64_t curr_vd = curr_wf_m.get_vd(idx);
            uint32_t node_id = internal::WavefrontArray::calc_node_id_from_vd(curr_vd);
            int32_t diag = internal::WavefrontArray::calc_diag_from_vd(curr_vd);
            int32_t j = curr_wf_m.get_offset(idx);
            int32_t i = diag + j;

            if (i + 1 == query_length &&
                j + 1 == static_cast<int32_t>(_patricia_tree.get_label_length(node_id))) {
                report_strings_at(node_id, current_score, results);
            }
        }

        // 停止条件チェック
        if (stop_predicate(current_score, results)) {
            break;
        }

        // Algorithm 3: DT-Patricia Expand
        expand(padded_query, wf_history_d, wf_history_m, wf_history_i, next_wf_array_d,
               next_wf_array_m, next_wf_array_i, curr_idx, history_size, _active_counts, buffer,
               pending_d, merged_wf_array_d, expand_scratch, _reached_d, current_score);
        uint32_t next_idx = internal::increment_mod(curr_idx, history_size);
        if (upper_bound >= 0) {
            prune_by_upper_bound<false>(
                wf_history_d[next_idx], wf_history_m[next_idx], wf_history_i[next_idx],
                _patricia_tree.get_subtree_max_lengths(), _patricia_tree.get_subtree_min_lengths(),
                query_length, upper_bound - (current_score + 1));
        }

        ++current_score;
        curr_idx = next_idx;
        next_idx = internal::increment_mod(next_idx, history_size);

        // 以降のループで使うことはないので履歴をリセット
        wf_history_d[next_idx].clear_logical_size();
        wf_history_m[next_idx].clear_logical_size();
        wf_history_i[next_idx].clear_logical_size();
    }
    state_guard.restore();
    return results;
}

}  // namespace dt_patricia
