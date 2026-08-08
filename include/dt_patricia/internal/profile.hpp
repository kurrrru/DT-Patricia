#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <ostream>

namespace dt_patricia::internal::profile {

// DT_PATRICIA_PROFILE が定義されているときだけ計数する。
// 未定義時は record_* が空になり、呼び出し側には何も残らない。
inline constexpr bool ENABLED =
#ifdef DT_PATRICIA_PROFILE
    true;
#else
    false;
#endif

// expand を密配列で行うかどうかの損益は、ノード 1 つあたりの
//   S = そのノードが持つ状態数
//   R = そのノードで埋める必要のある対角線範囲の幅
// の比で決まる。
//   - R/S が 1 に近ければ密配列が有利（fill がほぼ無駄にならない）
//   - S が小さいノードが多いほど、ノードあたりの固定費（fill と番兵ぶんの余白）が
//     相対的に重くなり、疎なマージのほうが速くなりうる
// 2 文字列の対角遷移と違い、木の上では波面が多数のノードに薄く分散するため、
// S の分布そのものが判断材料になる。
class Counters {
 public:
    static constexpr std::size_t STATE_BUCKET_COUNT = 9;
    static constexpr std::size_t SCORE_SLOT_COUNT = 64;

    void set_score(int32_t score) noexcept { _score = score; }

    void record_node(std::size_t states, std::size_t range) noexcept {
        ++_node_count;
        _state_total += states;
        _range_total += range;
        ++_state_hist[state_bucket(states)];

        const std::size_t slot = score_slot();
        ++_score_nodes[slot];
        _score_states[slot] += states;
        _score_range[slot] += range;
    }

    void record_extend_state(bool dominated) noexcept {
        ++_extend_states;
        if (dominated) {
            ++_extend_dominated;
        }
    }

    void reset() noexcept { *this = Counters(); }

    [[nodiscard]] std::size_t node_count() const noexcept { return _node_count; }
    [[nodiscard]] std::size_t state_total() const noexcept { return _state_total; }
    [[nodiscard]] std::size_t range_total() const noexcept { return _range_total; }

    void dump(std::ostream &os) const;

 private:
    // 1, 2, 3, 4, 5-8, 9-16, 17-32, 33-64, 65+
    [[nodiscard]] static std::size_t state_bucket(std::size_t states) noexcept {
        if (states <= 4) {
            return (states == 0) ? 0 : states - 1;
        }
        if (states <= 8) {
            return 4;
        }
        if (states <= 16) {
            return 5;
        }
        if (states <= 32) {
            return 6;
        }
        if (states <= 64) {
            return 7;
        }
        return 8;
    }

    [[nodiscard]] std::size_t score_slot() const noexcept {
        if (_score < 0) {
            return 0;
        }
        const std::size_t slot = static_cast<std::size_t>(_score);
        return (slot < SCORE_SLOT_COUNT) ? slot : SCORE_SLOT_COUNT - 1;
    }

    int32_t _score = 0;
    std::size_t _node_count = 0;
    std::size_t _state_total = 0;
    std::size_t _range_total = 0;
    std::array<std::size_t, STATE_BUCKET_COUNT> _state_hist = {};
    std::array<std::size_t, SCORE_SLOT_COUNT> _score_nodes = {};
    std::array<std::size_t, SCORE_SLOT_COUNT> _score_states = {};
    std::array<std::size_t, SCORE_SLOT_COUNT> _score_range = {};
    std::size_t _extend_states = 0;
    std::size_t _extend_dominated = 0;
};

inline void Counters::dump(std::ostream &os) const {
    os << "=== DT-Patricia expand profile ===\n";
    os << "nodes visited : " << _node_count << "\n";
    os << "states (sum S): " << _state_total << "\n";
    os << "range  (sum R): " << _range_total << "\n";
    if (_state_total > 0) {
        const double ratio = static_cast<double>(_range_total) / static_cast<double>(_state_total);
        os << "sum R / sum S : " << ratio << "\n";
    }
    if (_node_count > 0) {
        const double mean_s = static_cast<double>(_state_total) / static_cast<double>(_node_count);
        os << "mean S / node : " << mean_s << "\n";
    }

    if (_extend_states > 0) {
        const double rate =
            static_cast<double>(_extend_dominated) / static_cast<double>(_extend_states);
        os << "\nextend states  : " << _extend_states << "\n";
        os << "  dominated    : " << _extend_dominated << "  (" << (rate * 100.0) << "%)\n";
    }

    static constexpr std::array<const char *, STATE_BUCKET_COUNT> BUCKET_LABEL = {
        "1", "2", "3", "4", "5-8", "9-16", "17-32", "33-64", "65+"};
    os << "\nS histogram (nodes per bucket)\n";
    for (std::size_t i = 0; i < STATE_BUCKET_COUNT; ++i) {
        os << "  " << BUCKET_LABEL[i] << "\t" << _state_hist[i] << "\n";
    }

    os << "\nper score (slot " << (SCORE_SLOT_COUNT - 1) << " aggregates the tail)\n";
    os << "  score\tnodes\tsum S\tsum R\tR/S\n";
    for (std::size_t i = 0; i < SCORE_SLOT_COUNT; ++i) {
        if (_score_nodes[i] == 0) {
            continue;
        }
        os << "  " << i << "\t" << _score_nodes[i] << "\t" << _score_states[i] << "\t"
           << _score_range[i] << "\t";
        if (_score_states[i] > 0) {
            os << static_cast<double>(_score_range[i]) / static_cast<double>(_score_states[i]);
        }
        os << "\n";
    }
}

// 計測用のプロセス全体のカウンタ。
// [NOTE] 状態を引数で引き回さずに済ませるための意図的な例外。DT_PATRICIA_PROFILE 有効時
// 以外は更新されない。スレッド安全ではないので、プロファイルビルドは単一スレッドで回すこと。
inline Counters &counters() noexcept {
    static Counters instance;
    return instance;
}

inline void set_score(int32_t score) noexcept {
    if constexpr (ENABLED) {
        counters().set_score(score);
    }
}

inline void record_extend_state(bool dominated) noexcept {
    if constexpr (ENABLED) {
        counters().record_extend_state(dominated);
    }
}

inline void record_node(std::size_t states, std::size_t range) noexcept {
    if constexpr (ENABLED) {
        counters().record_node(states, range);
    }
}

}  // namespace dt_patricia::internal::profile
