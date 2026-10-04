#include <cassert>
#include <stdexcept>

#include <dt_patricia/aligner.hpp>

namespace dt_patricia {

// =========================================================
// 制限付き探索の設定と解除
// =========================================================

// 費用は O(N + U + U0)（U はこの制限で c が正になるノード数、U0 は前の制限のそれ）。
// 前の制限の後始末 (U0) は、解除ではなくここで払う。解除はクエリごとに行いうるので検索の時間に
// 含めるが、設定はクエリごとに行わなくて済むので含めない。
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::set_restriction(const std::vector<uint8_t> &mask) {
    const uint32_t string_count = _patricia_tree.string_count();
    if (mask.size() != string_count) {
        throw std::invalid_argument(
            "DTPatricia::set_restriction: mask size must equal string_count()");
    }
    discard_restriction();
    // 制限を使わない利用者に O(S) のメモリを持たせないよう、初回の設定時に確保する
    if (_restricted_counts.size() != _active_counts.size()) {
        _restricted_counts.assign(_active_counts.size(), 0);
    }

    try {
        _allowed = mask;
        const uint32_t root = _patricia_tree.root_id();
        uint32_t allowed = 0;
        for (uint32_t id = 0; id < string_count; ++id) {
            if (mask[id] == 0) {
                continue;
            }
            ++allowed;
            // 0 から 1 になったノードだけ親へ進む。各ノードが 0 から 1 になるのは一度だけ。
            uint32_t v = _patricia_tree.get_string_node(id);
            while (true) {
                if (_restricted_counts[v] != 0) {
                    ++_restricted_counts[v];
                    break;
                }
                // 先に積んでおけば、ここで例外が出ても discard_restriction() で戻せる
                _restricted_nodes.push_back(v);
                _restricted_counts[v] = 1;
                if (v == root) {
                    break;
                }
                v = _patricia_tree.get_parent(v);
            }
        }
        _restricted = true;
        _target_count = allowed;
    } catch (...) {
        discard_restriction();
        throw;
    }
#ifndef NDEBUG
    // 前の制限の値が残っても結果は変わらず枝刈りが弱くなるだけなので、結果の比較では見つからない。
    // ここで定義どおりの値と突き合わせる。
    assert(restricted_counts_match_definition());
#endif
}

// 制限の状態 (_restricted_counts, _allowed) には触れず、使わないようにするだけなので O(1)。
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::clear_restriction() noexcept {
    _restricted = false;
    _target_count = _patricia_tree.string_count();
}

// 直近の制限で書き換えた _restricted_counts を 0 に戻し、制限を解除する。O(U)。
// set_restriction からだけ呼ぶ（設定の費用として払う）。
template <AlphabetPolicy Alphabet, typename CostType>
void DTPatricia<Alphabet, CostType>::discard_restriction() noexcept {
    for (const uint32_t v : _restricted_nodes) {
        _restricted_counts[v] = 0;
    }
    _restricted_nodes.clear();
    clear_restriction();
}

// _restricted_counts が定義
//   c[v] = (v で終わる許可された ID の数) + (c[u] > 0 である終端以外の子 u の数)
// と一致するかを、ゼロから計算し直して調べる。O(S + N)。デバッグ用。
template <AlphabetPolicy Alphabet, typename CostType>
bool DTPatricia<Alphabet, CostType>::restricted_counts_match_definition() const {
    std::vector<uint32_t> expected(_restricted_counts.size(), 0);
    for (uint32_t id = 0; id < _patricia_tree.string_count(); ++id) {
        if (_allowed[id] != 0) {
            ++expected[_patricia_tree.get_string_node(id)];
        }
    }
    // ノード ID は BFS 順で親が子より小さいので、降順にたどれば子の値は親より先に確定する
    const uint32_t root = _patricia_tree.root_id();
    for (uint32_t u = static_cast<uint32_t>(expected.size()); u-- > root + 1;) {
        const uint32_t parent = _patricia_tree.get_parent(u);
        if (parent == 0 || _patricia_tree.transition(parent, tree_type::CODE_TERM) == u) {
            continue;  // 使われていない添字か、終端コードで入る子
        }
        if (expected[u] > 0) {
            ++expected[parent];
        }
    }
    return expected == _restricted_counts;
}

}  // namespace dt_patricia
