#include <cassert>

#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// =========================================================
// 探索の作業状態：枝刈り用の値 c の初期化、報告時の更新、探索後の復元
// =========================================================

// c を「v で終わる ID の数 + 終端以外の子の数」で埋める。
// Patricia 木では終端以外の子の部分木には必ず ID があるので、子の数がそのまま c[u] > 0 の数になる。
// 費用は O(S)（S は double array の添字空間の大きさ）。構築時に一度だけ行う。
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::init_active_counts() {
    const uint32_t node_count = _patricia_tree.node_count();
    _active_counts.assign(node_count, 0);
    for (uint32_t u = _patricia_tree.root_id() + 1; u < node_count; ++u) {
        const uint32_t parent = _patricia_tree.get_parent(u);
        if (parent == 0) {
            continue;  // 使われていない添字
        }
        if (_patricia_tree.transition(parent, tree_type::CODE_TERM) == u) {
            _active_counts[parent] += static_cast<uint32_t>(_patricia_tree.get_string_id(u).size());
        } else {
            _active_counts[parent] += 1;
        }
    }
}

// 費用は、報告した件数と、c が 0 になったノードの数に比例する。
// c が 0 になるノードは報告したノードの祖先なので到達ノードであり、1 回の探索で高々一度ずつ。
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::report_strings_at(uint32_t node_id, int32_t score,
                                                       std::vector<AlignmentResult> &results) {
    const uint32_t term_node = _patricia_tree.transition(node_id, tree_type::CODE_TERM);
    if (term_node == 0) {
        return;
    }
    uint32_t reported = 0;
    for (const uint32_t id : _patricia_tree.get_string_id(term_node)) {
        if (_reported[id] != 0 || (_restricted && _allowed[id] == 0)) {
            continue;
        }
        // results に積んでからフラグを立てる。逆順だと push_back が例外を投げたとき、
        // フラグだけ立って results に載らない ID ができ、復元で戻せなくなる。
        results.push_back({id, static_cast<uint32_t>(score)});
        _reported[id] = 1;
        ++reported;
    }
    if (reported == 0) {
        return;
    }

    std::vector<uint32_t> &counts = current_counts();
    const uint32_t root = _patricia_tree.root_id();
    _count_undo.emplace_back(node_id, counts[node_id]);
    assert(counts[node_id] >= reported);
    counts[node_id] -= reported;
    uint32_t v = node_id;
    while (counts[v] == 0 && v != root) {
        v = _patricia_tree.get_parent(v);
        _count_undo.emplace_back(v, counts[v]);
        assert(counts[v] > 0);
        counts[v] -= 1;
    }
}

// 費用は O(N_k + 書き換えた c の数 + 区画を確保したノードの数)。
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::restore_search_state(
    const std::vector<AlignmentResult> &results) noexcept {
    // _reported が 1 の ID は results に載った ID とちょうど一致する
    for (const AlignmentResult &r : results) {
        _reported[r.string_id] = 0;
    }
    // 同じノードを複数回書き換えた場合に最初の値が最後に書き戻されるよう、逆順にたどる
    std::vector<uint32_t> &counts = current_counts();
    for (auto it = _count_undo.rbegin(); it != _count_undo.rend(); ++it) {
        counts[it->first] = it->second;
    }
    _count_undo.clear();
    _reached.reset();
    _reached_d.reset();
}

}  // namespace dt_patricia
