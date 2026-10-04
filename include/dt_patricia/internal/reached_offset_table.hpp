#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace dt_patricia::internal {

// 支配判定をどこで行うかの切り替え。
//
//   既定（false）: 子ノードへプッシュする地点でのみ判定する
//   DT_PATRICIA_DOMINANCE_AT_ENTRY 定義時（true）: extend の入口で全状態を判定する
//
// 計算量のオーダー的にはDT_PATRICIA_DOMINANCE_AT_ENTRY 定義時の方が有利だが、実測ではfalse
// の方が速いのでfalseを既定とする
inline constexpr bool DOMINANCE_AT_EXTEND_ENTRY =
#ifdef DT_PATRICIA_DOMINANCE_AT_ENTRY
    true;
#else
    false;
#endif

// (ノード v, 対角線 d) ごとに「その状態が既に処理済みか」を記録する表。
//
// ラベル末尾でのみ判定する既定の設定では、判定地点のoffsetはlabel_len - 1であるため、
// その状態に到達しているかどうかの1ビットの記録で十分である。
// 一方、extend の入口で判定する設定では、状態ごとに異なる到達点を比べる必要があるため、
// 32 ビットの offset を保持しなければならない。
// この違いは、メンバ変数_arenaの各要素の使い方の違いによって実現する。
//
//   既定                             : 1 要素 = 対角線 32 本ぶんのビット
//   DT_PATRICIA_DOMINANCE_AT_ENTRY   : 1 要素 = 対角線 1 本ぶんの offset
class ReachedOffsetTable {
 public:
    // ノード数 0 の表。使わない層の表（linear の D 層）のためにある。
    // dominated() / record() を呼んではならない。
    ReachedOffsetTable() = default;

    // ノード添字で引く表 (_begin, _radius) をここで確保する。確保はノード数に比例するので、
    // クエリごとではなく所有者の構築時に一度だけ行い、クエリ間では reset() で使い回す。
    //   parent_path_len[v] : 根から v の親までのパス長。ウィンドウの中心になる。
    explicit ReachedOffsetTable(const std::vector<uint32_t> &parent_path_len)
        : _center(&parent_path_len),
          // 半径 0 を「このクエリではまだ確保していない」の印にする。
          // 実際の半径は必ず INITIAL_RADIUS 以上なので 0 と衝突しない。
          _begin(parent_path_len.size(), 0),
          _radius(parent_path_len.size(), 0) {}

    // クエリで書き換えた箇所だけを戻し、構築直後と同じ状態にする。
    // 費用はこのクエリで区画を確保したノードの数に比例する。
    void reset() noexcept {
        for (const uint32_t node_id : _allocated_nodes) {
            _radius[node_id] = 0;  // _begin は _radius が 0 のノードでは読まれないので戻さない
        }
        _allocated_nodes.clear();
        _arena.clear();  // 追加する要素は必ず ARENA_INIT で埋めるので、容量だけ残せばよい
        // 残すと、次のクエリの最初のノードが前回の最後と同じとき refresh_cache が省かれ、
        // 消した区画を指したまま使ってしまう
        _cached_node = INVALID_NODE;
    }

    // (node_id, diagonal) が記録済みの内容に照らして支配されているか。
    // 既定の設定では offset は使われず、記録の有無だけを返す。
    // score はウィンドウの不変条件 |d - L_p(v)| <= s の検査にのみ使う。
    [[nodiscard]] bool dominated(uint32_t node_id, int32_t diagonal, int32_t offset,
                                 int32_t score) {
        const std::size_t pos = locate(node_id, diagonal, score);
        if constexpr (DOMINANCE_AT_EXTEND_ENTRY) {
            return offset <= static_cast<int32_t>(_arena[_cached_begin + pos]);
        } else {
            const uint32_t word = _arena[_cached_begin + pos / BITS_PER_WORD];
            return ((word >> (pos % BITS_PER_WORD)) & 1u) != 0u;
        }
    }

    // (node_id, diagonal) を処理済みとして記録する。
    // 既定の設定では offset を無視してビットを立てるだけ。
    void record(uint32_t node_id, int32_t diagonal, int32_t offset, int32_t score) {
        const std::size_t pos = locate(node_id, diagonal, score);
        if constexpr (DOMINANCE_AT_EXTEND_ENTRY) {
            _arena[_cached_begin + pos] = static_cast<uint32_t>(offset);
        } else {
            _arena[_cached_begin + pos / BITS_PER_WORD] |= (uint32_t{1} << (pos % BITS_PER_WORD));
        }
    }

 private:
    static constexpr std::size_t BITS_PER_WORD = 32;
    static constexpr uint32_t INVALID_NODE = UINT32_MAX;

    // ビット詰めの伸長でワード単位のコピーに収めるため、BITS_PER_WORD の倍数にする。
    static constexpr std::size_t INITIAL_RADIUS = BITS_PER_WORD;

    // 未到達を表す値。
    static constexpr uint32_t ARENA_INIT =
        DOMINANCE_AT_EXTEND_ENTRY ? static_cast<uint32_t>(INT32_MIN) : 0u;

    // 半径 radius（対角線の本数）の区画が占めるアリーナの要素数。
    // radius は BITS_PER_WORD の倍数なので、ビット詰めでも割り切れる。
    [[nodiscard]] static constexpr std::size_t elements_for(std::size_t radius) noexcept {
        if constexpr (DOMINANCE_AT_EXTEND_ENTRY) {
            return radius * 2;
        } else {
            return (radius * 2) / BITS_PER_WORD;
        }
    }

    // 対角線 shift 本ぶんのずれが、アリーナの要素いくつぶんにあたるか。
    [[nodiscard]] static constexpr std::size_t elements_shift(std::size_t shift) noexcept {
        if constexpr (DOMINANCE_AT_EXTEND_ENTRY) {
            return shift;
        } else {
            return shift / BITS_PER_WORD;
        }
    }

    // 区画内での対角線の通し番号を返す。未確保なら確保し、範囲外なら区画を伸ばす。
    // 呼び出し後は _cached_begin がそのノードの区画先頭を指している。
    [[nodiscard]] std::size_t locate(uint32_t node_id, int32_t diagonal,
                                     [[maybe_unused]] int32_t score) {
        // extend は状態をノードごとにまとめて処理するので、直前のノードの区画情報を
        // 持ち回れば、状態 1 個あたりの表引きが 4 本から 1 本に減る。
        if (node_id != _cached_node) {
            refresh_cache(node_id);
        }

        const int64_t delta = static_cast<int64_t>(diagonal) - _cached_center;
        assert(delta >= -static_cast<int64_t>(score) && delta <= static_cast<int64_t>(score) &&
               "|d - L_p(v)| <= s に違反");

        if (delta < -_cached_radius || delta >= _cached_radius) [[unlikely]] {
            grow(node_id, delta);
            refresh_cache(node_id);
        }

        // 区画は [center - radius, center + radius) を覆うので、
        // delta に radius を足すと区画の先頭から数えた通し番号になる。
        return static_cast<std::size_t>(delta + _cached_radius);
    }

    // ノードの区画情報をキャッシュに読み込む。未確保なら確保する。
    // アリーナは伸長で再確保されうるので、ポインタではなく添字を持つ。
    void refresh_cache(uint32_t node_id) {
        if (_radius[node_id] == 0) {
            allocate(node_id, 0);
        }
        _cached_node = node_id;
        _cached_center = static_cast<int64_t>((*_center)[node_id]);
        _cached_begin = _begin[node_id];
        _cached_radius = static_cast<int64_t>(_radius[node_id]);
    }

    [[nodiscard]] static std::size_t radius_covering(int64_t delta) noexcept {
        const int64_t need = (delta < 0) ? -delta : delta;
        std::size_t radius = INITIAL_RADIUS;
        while (static_cast<int64_t>(radius) <= need) {
            radius *= 2;
        }
        return radius;
    }

    void allocate(uint32_t node_id, int64_t delta) {
        const std::size_t radius = radius_covering(delta);
        // 先に積んでおけば、後続の確保が例外を投げても reset() で戻せる
        _allocated_nodes.push_back(node_id);
        _begin[node_id] = _arena.size();
        _radius[node_id] = static_cast<uint32_t>(radius);
        _arena.resize(_arena.size() + elements_for(radius), ARENA_INIT);
    }

    void grow(uint32_t node_id, int64_t delta) {
        const std::size_t old_radius = _radius[node_id];
        const std::size_t old_begin = _begin[node_id];
        const std::size_t new_radius = radius_covering(delta);

        const std::size_t new_begin = _arena.size();
        _arena.resize(_arena.size() + elements_for(new_radius), ARENA_INIT);
        const std::size_t shift = elements_shift(new_radius - old_radius);
        const std::size_t old_elements = elements_for(old_radius);
        for (std::size_t i = 0; i < old_elements; ++i) {
            _arena[new_begin + shift + i] = _arena[old_begin + i];
        }
        _begin[node_id] = new_begin;
        _radius[node_id] = static_cast<uint32_t>(new_radius);
    }

    const std::vector<uint32_t> *_center = nullptr;
    std::vector<uint32_t> _arena;
    std::vector<std::size_t> _begin;  // _arenaにおける区間の開始位置。インデックスはnode_id
    std::vector<uint32_t> _radius;  // 確保済み領域の半径。インデックスはnode_id
    std::vector<uint32_t> _allocated_nodes;  // このクエリで区画を確保したノード。reset() で使う

    // キャッシュ情報
    uint32_t _cached_node = INVALID_NODE;
    std::size_t _cached_begin = 0;
    int64_t _cached_center = 0;
    int64_t _cached_radius = 0;
};

}  // namespace dt_patricia::internal
