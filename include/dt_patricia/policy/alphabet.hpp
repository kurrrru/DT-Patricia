#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>

namespace dt_patricia {

namespace detail {

// alphabet policy が満たすべき不変条件をコンパイル時に検証する。
// - CODE_TERM は必ず 0
// - CODE_MAX は 0 より大きい
// - make_char_to_code()['\0'] == CODE_TERM
// - make_char_to_code() の全要素は 0..CODE_MAX に収まる
template <class Alphabet>
consteval bool validate_alphabet_policy() {
    if (Alphabet::CODE_TERM != 0) {
        return false;
    }
    if (Alphabet::CODE_MAX == 0) {
        return false;
    }

    constexpr auto table = Alphabet::make_char_to_code();

    if (table['\0'] != Alphabet::CODE_TERM) {
        return false;
    }
    for (uint8_t code : table) {
        if (code > Alphabet::CODE_MAX) {
            return false;
        }
    }
    return true;
}

}  // namespace detail

// alphabet policy の公開要件。
// policy は CHAR_TO_CODE 実体を持たず、consteval table generator を提供するだけにする。
// CODE_A / CODE_C / ... のような個別コードは公開要件に含めない（各 policy の
// make_char_to_code() 内の local constexpr に閉じ込める）。
template <class Alphabet>
concept AlphabetPolicy =
    requires {
        { Alphabet::CODE_TERM } -> std::convertible_to<uint8_t>;
        { Alphabet::CODE_MAX } -> std::convertible_to<uint8_t>;
        { Alphabet::make_char_to_code() } noexcept
            -> std::same_as<std::array<uint8_t, 256>>;
    } &&
    detail::validate_alphabet_policy<Alphabet>();

// =========================================================
// DNA alphabet: A/C/G/T (U は T として扱う) + その他
// =========================================================
struct DnaAlphabet {
    static constexpr uint8_t CODE_TERM = 0;
    static constexpr uint8_t CODE_MAX  = 5;

    static consteval std::array<uint8_t, 256> make_char_to_code() noexcept {
        constexpr uint8_t CODE_A     = 1;
        constexpr uint8_t CODE_C     = 2;
        constexpr uint8_t CODE_G     = 3;
        constexpr uint8_t CODE_T     = 4;
        constexpr uint8_t CODE_OTHER = 5;

        std::array<uint8_t, 256> table{};
        table.fill(CODE_OTHER);

        table['\0'] = CODE_TERM;

        table['A'] = CODE_A; table['a'] = CODE_A;
        table['C'] = CODE_C; table['c'] = CODE_C;
        table['G'] = CODE_G; table['g'] = CODE_G;
        table['T'] = CODE_T; table['t'] = CODE_T;

        // RNA の U は T として扱う
        table['U'] = CODE_T; table['u'] = CODE_T;

        // N/n やその他の文字はすべて table.fill(CODE_OTHER) により OTHER になる

        return table;
    }
};

// =========================================================
// R/Y alphabet: purine(R) / pyrimidine(Y) + その他
// =========================================================
// R/Y 表記済み文字列にも、DNA 文字列を R/Y に潰す用途にも対応する。
struct RyAlphabet {
    static constexpr uint8_t CODE_TERM = 0;
    static constexpr uint8_t CODE_MAX  = 3;

    static consteval std::array<uint8_t, 256> make_char_to_code() noexcept {
        constexpr uint8_t CODE_R     = 1;
        constexpr uint8_t CODE_Y     = 2;
        constexpr uint8_t CODE_OTHER = 3;

        std::array<uint8_t, 256> table{};
        table.fill(CODE_OTHER);

        table['\0'] = CODE_TERM;

        // purine
        table['A'] = CODE_R; table['a'] = CODE_R;
        table['G'] = CODE_R; table['g'] = CODE_R;
        table['R'] = CODE_R; table['r'] = CODE_R;

        // pyrimidine
        table['C'] = CODE_Y; table['c'] = CODE_Y;
        table['T'] = CODE_Y; table['t'] = CODE_Y;
        table['U'] = CODE_Y; table['u'] = CODE_Y;
        table['Y'] = CODE_Y; table['y'] = CODE_Y;

        // unknown / N / その他文字は CODE_OTHER に落ちる

        return table;
    }
};

static_assert(AlphabetPolicy<DnaAlphabet>);
static_assert(AlphabetPolicy<RyAlphabet>);

}  // namespace dt_patricia
