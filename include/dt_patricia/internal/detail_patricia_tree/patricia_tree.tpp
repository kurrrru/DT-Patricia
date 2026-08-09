#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <dt_patricia/patricia_tree.hpp>

namespace dt_patricia {

namespace detail_patricia_tree {
// MSD Radix Sortの実装
// depth: 現在見ている文字位置
// [start, end): ソート対象のインデックス範囲
template <AlphabetPolicy Alphabet>
void msd_radix_sort_recursive(std::vector<uint32_t> &indices,
                              const std::vector<std::string> &input_data, size_t start, size_t end,
                              size_t depth);
}  // namespace detail_patricia_tree

template <AlphabetPolicy Alphabet>
PatriciaTree<Alphabet>::PatriciaTree(const std::vector<std::string> &input_data)
    : _size(static_cast<uint32_t>(input_data.size())) {
    if (input_data.empty()) {
        return;
    }

    // MSD Radix Sortでインデックスをソート
    std::vector<uint32_t> indices(input_data.size());
    std::iota(indices.begin(), indices.end(), 0);
    detail_patricia_tree::msd_radix_sort_recursive<Alphabet>(indices, input_data, 0, indices.size(),
                                                             0);

    build(input_data, indices);
}

template <AlphabetPolicy Alphabet>
void PatriciaTree<Alphabet>::build(const std::vector<std::string> &input_data,
                                   const std::vector<uint32_t> &sorted_indices) {
    // ===========================================================
    // 1. ツリー構築のための初期化
    // ===========================================================
    _base.reserve(input_data.size() * 2);
    _check.reserve(input_data.size() * 2);
    _subtree_counts.reserve(input_data.size() * 2);

    // ノードID 0 は無効値、ノードID 1 をルートノードとする
    _base.assign(2, 0);
    _check.assign(2, 0);
    _label_offset.assign(2, 0);
    _label_len.assign(2, 0);
    _string_ids_offset.assign(2, 0);
    _string_ids_count.assign(2, 0);
    _subtree_counts.assign(2, 0);

    if (sorted_indices.size() > 0) {
        _subtree_counts[1] = static_cast<uint32_t>(sorted_indices.size());
    }

    // ===========================================================
    // 2. BFSでツリー構築
    // ===========================================================
    // BFS用のキュー
    std::vector<QueueItem> queue;
    queue.reserve(input_data.size());
    queue.push_back({1, 0, sorted_indices.size(), 0});

    size_t queue_pos = 0;

    while (queue_pos < queue.size()) {
        const auto item = queue[queue_pos++];
        build_node_bfs(item.node_id, item.start_idx, item.end_idx, item.label_offset, input_data,
                       sorted_indices, queue);
    }
    // SIMDのためのパディング
    TEXT_POOL.append(SIMD_PADDING_SIZE, '\0');

    // ===========================================================
    // 3. 部分木の長さの最小/最大の事前計算
    // ===========================================================
    compute_subtree_length_bounds();
    compute_parent_path_lengths();
}

// 根から各ノードの親までのパス長。ノード ID は BFS 順で親が必ず子より小さいので、
// 昇順の 1 パスで確定する。
template <AlphabetPolicy Alphabet>
void PatriciaTree<Alphabet>::compute_parent_path_lengths() {
    const uint32_t num_nodes = static_cast<uint32_t>(_base.size());
    _parent_path_len.assign(num_nodes, 0);
    std::vector<uint32_t> path_len(num_nodes, 0);

    const uint32_t root = root_id();
    if (num_nodes <= root) {
        return;
    }
    path_len[root] = _label_len[root];
    for (uint32_t v = root + 1; v < num_nodes; ++v) {
        const uint32_t parent = _check[v];
        if (parent == 0) {
            continue;
        }
        _parent_path_len[v] = path_len[parent];
        path_len[v] = path_len[parent] + _label_len[v];
    }
}

template <AlphabetPolicy Alphabet>
void PatriciaTree<Alphabet>::compute_subtree_length_bounds() {
    uint32_t num_nodes = static_cast<uint32_t>(_base.size());
    _subtree_min_len.assign(num_nodes, 0xFFFFFFFF);
    _subtree_max_len.assign(num_nodes, 0);

    // 2. 逆順ループによるボトムアップ集計（第1パスのみ）
    // ID 1がルートなので、1まで処理する。0は無効値なので飛ばす。
    for (uint32_t i = num_nodes - 1; i >= 1; --i) {
        if (_check[i] == 0 && i != root_id()) {
            continue;
        }
        uint32_t label_l = _label_len[i];

        if (is_leaf(i)) {
            _subtree_min_len[i] = label_l;
            _subtree_max_len[i] = label_l;
        } else {
            _subtree_max_len[i] += label_l;
            _subtree_min_len[i] += label_l;
        }

        // 親ノードへの伝播
        uint32_t p = _check[i];
        if (p != 0) {
            uint32_t child_full_min = _subtree_min_len[i];
            uint32_t child_full_max = _subtree_max_len[i];

            if (child_full_min < _subtree_min_len[p]) {
                _subtree_min_len[p] = child_full_min;
            }
            if (child_full_max > _subtree_max_len[p]) {
                _subtree_max_len[p] = child_full_max;
            }
        }
    }
}

template <AlphabetPolicy Alphabet>
void PatriciaTree<Alphabet>::build_node_bfs(uint32_t node_id, size_t start_idx, size_t end_idx,
                                            size_t label_offset,
                                            const std::vector<std::string> &input_data,
                                            const std::vector<uint32_t> &sorted_indices,
                                            std::vector<QueueItem> &queue) {
    // 1. ラベル（共通接頭辞）の決定
    const std::string &first_str = input_data[sorted_indices[start_idx]];
    const std::string &last_str = input_data[sorted_indices[end_idx - 1]];
    size_t max_lcp = std::min(first_str.size(), last_str.size()) - label_offset;
    size_t lcp_len = 0;
    while (lcp_len < max_lcp &&
           CHAR_TO_CODE[static_cast<unsigned char>(first_str[label_offset + lcp_len])] ==
               CHAR_TO_CODE[static_cast<unsigned char>(last_str[label_offset + lcp_len])]) {
        ++lcp_len;
    }

    _label_offset[node_id] = static_cast<uint32_t>(TEXT_POOL.size());
    _label_len[node_id] = static_cast<uint32_t>(lcp_len);

    if (lcp_len > 0) {
        size_t pool_pos = TEXT_POOL.size();
        TEXT_POOL.append(first_str, label_offset, lcp_len);
        canonicalize_inplace<Alphabet>(TEXT_POOL.data() + pool_pos, lcp_len);
    }

    size_t current_offset = label_offset + lcp_len;

    // 2. グループ化
    struct Group {
        uint8_t code;
        size_t start;
        size_t end;
    };
    std::vector<Group> groups;

    auto get_code = [&](size_t idx) -> uint8_t {
        const std::string &s = input_data[sorted_indices[idx]];
        if (current_offset >= s.size()) {
            return CODE_TERM;
        }
        return CHAR_TO_CODE[static_cast<unsigned char>(s[current_offset])];
    };

    size_t group_start = start_idx;
    uint8_t current_code = get_code(group_start);

    for (size_t i = start_idx + 1; i < end_idx; ++i) {
        uint8_t code = get_code(i);
        if (code != current_code) {
            groups.push_back({current_code, group_start, i});
            current_code = code;
            group_start = i;
        }
    }
    groups.push_back({current_code, group_start, end_idx});

    // 3. ベース値の決定と配列リサイズ
    uint32_t base = static_cast<uint32_t>(_check.size());
    _base[node_id] = base;

    uint8_t max_code = 0;
    for (const auto &g : groups) {
        if (g.code > max_code) {
            max_code = g.code;
        }
    }
    size_t required_size = base + max_code + 1;

    if (_check.size() < required_size) {
        _base.resize(required_size, 0);
        _check.resize(required_size, 0);
        _label_offset.resize(required_size, 0);
        _label_len.resize(required_size, 0);
        _string_ids_offset.resize(required_size, 0);
        _string_ids_count.resize(required_size, 0);
        _subtree_counts.resize(required_size, 0);
    }

    // 4. 子ノードの作成（BFS: キューに追加）
    for (const auto &group : groups) {
        uint32_t child_node_id = base + group.code;
        _check[child_node_id] = node_id;

        // 部分木の単語数を更新
        uint32_t group_size = static_cast<uint32_t>(group.end - group.start);
        _subtree_counts[child_node_id] = group_size;

        if (group.code == CODE_TERM) {
            // 葉ノード
            _base[child_node_id] = 0;

            uint32_t ids_start = static_cast<uint32_t>(_string_ids.size());
            uint32_t ids_count = static_cast<uint32_t>(group.end - group.start);

            _string_ids_offset[child_node_id] = ids_start;
            _string_ids_count[child_node_id] = ids_count;

            _string_ids.insert(_string_ids.end(), sorted_indices.begin() + group.start,
                               sorted_indices.begin() + group.end);

            _label_offset[child_node_id] = 0;
            _label_len[child_node_id] = 0;
        } else {
            // 内部ノード: キューに追加
            queue.push_back({child_node_id, group.start, group.end, current_offset});
        }
    }
}

namespace detail_patricia_tree {
// MSD Radix Sort(再帰と非再帰のハイブリッド)
// - 最大バケットは再帰せず while ループで反復 → 再帰深さは文字列長に依存せず O(log n)
// - temp はスコープを閉じて降下前に解放 → 同時ヒープを O(n) に抑制
// - 共通接頭辞（単一バケット）は分配/コピーを省いて depth を進める
// depth: 現在見ている文字位置 / [start, end): ソート対象のインデックス範囲
template <AlphabetPolicy Alphabet>
void msd_radix_sort_recursive(std::vector<uint32_t> &indices,
                              const std::vector<std::string> &input_data, size_t start, size_t end,
                              size_t depth) {
    using Tree = PatriciaTree<Alphabet>;

    constexpr size_t INSERTION_SORT_THRESHOLD = 16;
    constexpr size_t BUCKET_SIZE = Tree::BUCKET_SIZE;

    // 最大バケットはこの while ループで [start, end, depth] を更新して反復処理し、
    // 再帰するのは「最大以外」のバケットだけにする。
    while (true) {
        // --- 小区間は比較ソートに委譲 ---
        if (end - start <= INSERTION_SORT_THRESHOLD) {
            std::sort(indices.begin() + start, indices.begin() + end, [&](uint32_t a, uint32_t b) {
                const std::string &sa = input_data[a];
                const std::string &sb = input_data[b];
                size_t pos = depth;
                // コードで比較（MSD Radix Sortと一貫性を保つ）
                while (pos < sa.size() && pos < sb.size()) {
                    uint8_t code_a = Tree::CHAR_TO_CODE[static_cast<unsigned char>(sa[pos])];
                    uint8_t code_b = Tree::CHAR_TO_CODE[static_cast<unsigned char>(sb[pos])];
                    if (code_a != code_b) {
                        return code_a < code_b;
                    }
                    ++pos;
                }
                // 一方が他方のプレフィックスなら短い方（終端）が先
                return sa.size() < sb.size();
            });
            return;
        }

        // --- 1. カウント（この時点で count[c+1] = コード c の個数） ---
        std::array<size_t, BUCKET_SIZE + 1> count = {};
        for (size_t i = start; i < end; ++i) {
            const std::string &s = input_data[indices[i]];
            uint8_t code = (depth >= s.size())
                               ? Tree::CODE_TERM
                               : Tree::CHAR_TO_CODE[static_cast<unsigned char>(s[depth])];
            count[code + 1]++;
        }

        // --- 2. 単一バケット（共通接頭辞スキップ） ---
        // 非空バケットが1つだけなら分配もコピーも不要。
        //   終端のみ  → 全要素が等しい（ここで終わる）ので完了。
        //   非終端のみ → 全員が同じ文字で続く=共通接頭辞。depth を進めて反復（再帰しない）。
        {
            size_t nonempty = 0;
            size_t sole_code = 0;
            for (size_t c = 0; c < BUCKET_SIZE; ++c) {
                if (count[c + 1] > 0) {
                    ++nonempty;
                    sole_code = c;
                }
            }
            if (nonempty == 1) {
                if (sole_code == Tree::CODE_TERM) {
                    return;
                }
                ++depth;
                continue;
            }
        }

        // --- 3. 累積和（各バケットの開始位置） ---
        for (size_t i = 0; i < BUCKET_SIZE; ++i) {
            count[i + 1] += count[i];
        }
        const std::array<size_t, BUCKET_SIZE + 1> bucket_bounds = count;

        // --- 4. 分配とコピーバック（temp はここで確保し、降下前に解放する） ---
        {
            std::vector<uint32_t> temp(end - start);
            for (size_t i = start; i < end; ++i) {
                const std::string &s = input_data[indices[i]];
                uint8_t code = (depth >= s.size())
                                   ? Tree::CODE_TERM
                                   : Tree::CHAR_TO_CODE[static_cast<unsigned char>(s[depth])];
                temp[count[code]++] = indices[i];
            }
            std::copy(temp.begin(), temp.end(), indices.begin() + start);
        }

        // --- 5. 最大バケットは反復、それ以外は再帰 ---
        // 終端バケットは全要素が等しいので処理不要。非終端かつ2要素以上のバケットのみ対象。
        size_t best_start = 0;
        size_t best_end = 0;
        size_t best_size = 0;
        for (size_t c = 0; c < BUCKET_SIZE; ++c) {
            if (c == Tree::CODE_TERM) {
                continue;
            }
            size_t b_start = start + bucket_bounds[c];
            size_t b_end = start + bucket_bounds[c + 1];
            if (b_end - b_start < 2) {
                continue;  // 0/1 要素は既にソート済み
            }
            if (b_end - b_start > best_size) {
                // それまでの最大候補は再帰に回し、新しい最大を採用
                if (best_size >= 2) {
                    msd_radix_sort_recursive<Alphabet>(indices, input_data, best_start, best_end,
                                                       depth + 1);
                }
                best_start = b_start;
                best_end = b_end;
                best_size = b_end - b_start;
            } else {
                msd_radix_sort_recursive<Alphabet>(indices, input_data, b_start, b_end, depth + 1);
            }
        }

        // 最大バケットだけは反復で処理（再帰しない）
        if (best_size >= 2) {
            start = best_start;
            end = best_end;
            ++depth;
            continue;
        }
        return;
    }
}

}  // namespace detail_patricia_tree

}  // namespace dt_patricia
