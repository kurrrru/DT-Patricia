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
#include <utility>
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

    // 探索で使う作業領域（ノード数・文字列数に比例する）はここで一度だけ確保する。
    // 作業領域は探索中に書き換わるので、並列に検索する場合はスレッドごとに DTPatricia を
    // 一つ持つこと。索引 (PatriciaTree) はスレッド間で共有してよい。
    DTPatricia(const tree_type &patricia_tree, CostType cost = CostType())
        : _patricia_tree(patricia_tree),
          _cost(cost),
          _target_count(patricia_tree.string_count()),
          _reported(patricia_tree.string_count(), 0),
          _reached(patricia_tree.get_parent_path_lengths()),
          _reached_d(CostType::is_linear
                         ? internal::ReachedOffsetTable()
                         : internal::ReachedOffsetTable(patricia_tree.get_parent_path_lengths())) {
        init_active_counts();
    }
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
    // 3. 制限付き探索の設定
    // =========================================================

    // 以降の探索の対象を、mask[id] != 0 である文字列 ID に限る。
    // mask の長さは string_count() と等しくなければならない（違えば std::invalid_argument）。
    // 費用は O(N)。初回だけノード数に比例する確保が加わる。
    void set_restriction(const std::vector<uint8_t> &mask);

    // 制限を解除し、全 ID を探索対象に戻す。費用は O(1)。
    // 制限で書き換えた箇所は残すが、制限がなければ読まれない。次の set_restriction で消す。
    void clear_restriction() noexcept;

    [[nodiscard]] bool restricted() const noexcept { return _restricted; }

    // =========================================================
    // 4. アラインメントAPI
    // =========================================================
    // 制限が設定されていれば、その対象の中だけを探索する。

    std::vector<AlignmentResult> ed_to_all(const std::string &query) {
        std::size_t max_results = _target_count;
        return search_kernel(
            query, [max_results](int, const auto &r) { return r.size() == max_results; }, -1);
    }

    std::vector<AlignmentResult> ed_within_k(const std::string &query, int k) {
        std::size_t max_results = _target_count;
        return search_kernel(
            query,
            [k, max_results](int current_score, const auto &r) {
                return current_score >= k || r.size() == max_results;
            },
            k);
    }

    std::vector<AlignmentResult> ed_pth_smallest(const std::string &query, size_t p) {
        if (p > _target_count) {
            p = _target_count;
        }
        return search_kernel(query, [p](int, const auto &r) { return r.size() >= p; }, -1);
    }

    template <typename StopPredicate>
    std::vector<AlignmentResult> search_kernel(const std::string &query,
                                               StopPredicate stop_predicate, int upper_bound = -1)
        requires(CostType::is_linear);

    template <typename StopPredicate>
    std::vector<AlignmentResult> search_kernel(const std::string &query,
                                               StopPredicate stop_predicate, int upper_bound = -1)
        requires(!CostType::is_linear);

 private:
    const tree_type &_patricia_tree;
    CostType _cost;

    // 枝刈り用の値 c[v]。
    //   c[v] = (v で終わる ID のうち、探索対象で未報告のものの数) + (c[u] > 0 である子 u の数)
    // 子には終端コードで入る子を含めない。c[v] > 0 と「v の部分木に探索対象で未報告の ID がある」
    // は同値なので、枝刈りは 0 かどうかだけを見る。探索の外では構築直後（制限付きなら設定直後）の
    // 値を保つ。
    std::vector<uint32_t> _active_counts;  // 制限なしのときの c
    // 制限付きのときの c。_restricted が false のときは読まれない。
    // _restricted_nodes に載っていないノードは常に 0。
    std::vector<uint32_t> _restricted_counts;
    std::vector<uint32_t> _restricted_nodes;  // 直近の制限で _restricted_counts を書き換えたノード
    bool _restricted = false;
    uint32_t _target_count;  // 探索対象の ID の数

    // ID ごとに、この探索で報告済みなら 1。探索の外では全要素 0。
    std::vector<uint8_t> _reported;
    // ID ごとに、直近の制限で許可されていれば 0 以外。_restricted が false のときは読まれない。
    // c だけでは、許可した ID の祖先にあたるノードで終わる、許可していない ID を除けないので、
    // 報告時にこれで絞る。
    std::vector<uint8_t> _allowed;

    // 探索中に書き換えた c の（ノード, 書き換え前の値）。探索後に逆順に書き戻す。
    std::vector<std::pair<uint32_t, uint32_t>> _count_undo;

    internal::ReachedOffsetTable _reached;  // M 層（linear では唯一）の子生成の重複を潰す表
    internal::ReachedOffsetTable _reached_d;  // affine の D 層の子生成の重複を潰す表

    // 探索で書き換えた状態を、探索前に戻す。例外で探索を抜けた場合もデストラクタで戻す。
    class SearchStateGuard {
     public:
        SearchStateGuard(DTPatricia &aligner, const std::vector<AlignmentResult> &results) noexcept
            : _aligner(aligner), _results(results) {}
        ~SearchStateGuard() { restore(); }
        SearchStateGuard(const SearchStateGuard &) = delete;
        SearchStateGuard &operator=(const SearchStateGuard &) = delete;

        // 正常終了時は results を返す（ムーブする）前にこれを呼ぶ。
        // デストラクタだけに任せると、NRVO が効かない場合に空の results を見てしまう。
        void restore() noexcept {
            if (!_restored) {
                _aligner.restore_search_state(_results);
                _restored = true;
            }
        }

     private:
        DTPatricia &_aligner;
        const std::vector<AlignmentResult> &_results;
        bool _restored = false;
    };

    void init_active_counts();
    void restore_search_state(const std::vector<AlignmentResult> &results) noexcept;
    void discard_restriction() noexcept;
    [[nodiscard]] bool restricted_counts_match_definition() const;

    // 探索が読み書きする c。制限があれば制限付きのもの。
    [[nodiscard]] std::vector<uint32_t> &current_counts() noexcept {
        return _restricted ? _restricted_counts : _active_counts;
    }

    // ノード node_id で終わる ID のうち、探索対象で未報告のものを score で報告し、c を更新する。
    void report_strings_at(uint32_t node_id, int32_t score, std::vector<AlignmentResult> &results);

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
#include <dt_patricia/internal/detail_aligner/restriction.tpp>
#include <dt_patricia/internal/detail_aligner/search_kernel.tpp>
#include <dt_patricia/internal/detail_aligner/search_state.tpp>
