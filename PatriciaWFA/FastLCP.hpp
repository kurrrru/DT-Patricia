#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <bit>
#include <algorithm>

// ハードウェア固有のSIMDヘッダのインクルード
#if defined(__AVX2__)
    #include <immintrin.h>
#elif defined(__aarch64__) || defined(__ARM_NEON) || defined(_M_ARM64)
    #include <arm_neon.h>
#endif

/**
 * @brief パディングを前提とした完全無分岐・SIMD最適化LCP (Longest Common Prefix) 計算
 * @param s1 比較文字列1のポインタ
 * @param s2 比較文字列2のポインタ
 * @param max_len 比較する最大長（これを超えた一致は切り捨てられる）
 * @return 一致したバイト数 (最大 max_len)
 * @warning 文字列の末尾には少なくとも32バイトのパディング（ダミーデータ）が
 * 割り当てられている領域が続くことを前提とします。
 */
inline size_t fast_lcp(const char* s1, const char* s2, size_t max_len) {
    size_t matched = 0;

    if (max_len < 8) {
        while (matched < max_len && s1[matched] == s2[matched]) {
            matched++;
        }
        return matched;
    }

#if defined(__AVX2__)
    // =========================================================
    // Windows / x86_64 環境 (AVX2: 32バイト一括処理)
    // =========================================================
    while (matched < max_len) {
        // パディングがあるため、max_lenの残りが32バイト未満でも安全にロード可能
        __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(s1 + matched));
        __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(s2 + matched));
        
        __m256i cmp = _mm256_cmpeq_epi8(v1, v2);
        
        // _mm256_movemask_epi8: 各バイトの最上位ビットを集めて32bitのマスクを作成
        // 一致していれば0xFFになるため、マスクのビットは1になる。
        // ビット反転(~)により、最初に不一致となったバイトの位置に1が立つ。
        uint32_t mask = ~static_cast<uint32_t>(_mm256_movemask_epi8(cmp));
        
        if (mask != 0) {
            // 右から連続する0の数 ＝ 最初に1が立ったビット位置 ＝ 一致したバイト数
            matched += std::countr_zero(mask);
            return std::min(matched, max_len);
        }
        matched += 32;
    }

#elif defined(__aarch64__) || defined(__ARM_NEON) || defined(_M_ARM64)
    // =========================================================
    // Mac / ARM64 環境 (NEON: 16バイト一括処理)
    // =========================================================
    while (matched < max_len) {
        uint8x16_t v1 = vld1q_u8(reinterpret_cast<const uint8_t*>(s1 + matched));
        uint8x16_t v2 = vld1q_u8(reinterpret_cast<const uint8_t*>(s2 + matched));
        
        // vceqq_u8: 各バイトごとに一致で 0xFF, 不一致で 0x00
        uint8x16_t cmp = vceqq_u8(v1, v2);
        
        // NEONには x86 の movemask に直接該当する命令がないため、
        // 64ビット整数2つとしてキャストし、レーンごとに検証する。
        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        
        uint64_t lane0 = vgetq_lane_u64(cmp64, 0);
        // ~0ULL (すべてのビットが1) でなければ、どこかに不一致(0x00)が含まれている
        if (lane0 != ~0ULL) {
            // リトルエンディアン前提: 不一致バイト(0xFFに反転)の最初の位置を8で割りバイト数へ
            matched += std::countr_zero(~lane0) / 8;
            return std::min(matched, max_len);
        }
        
        uint64_t lane1 = vgetq_lane_u64(cmp64, 1);
        if (lane1 != ~0ULL) {
            matched += 8 + (std::countr_zero(~lane1) / 8);
            return std::min(matched, max_len);
        }
        
        matched += 16;
    }

#else
    // =========================================================
    // フォールバック (SWAR: 8バイト一括処理)
    // AVX2もNEONも有効でない場合の安全網
    // =========================================================
    while (matched < max_len) {
        uint64_t v1, v2;
        std::memcpy(&v1, s1 + matched, sizeof(uint64_t));
        std::memcpy(&v2, s2 + matched, sizeof(uint64_t));
        
        uint64_t diff = v1 ^ v2;
        if (diff != 0) {
            // XORで差分を取ると、一致したビットは0、異なるビットは1になる
            matched += std::countr_zero(diff) / 8;
            return std::min(matched, max_len);
        }
        matched += sizeof(uint64_t);
    }
#endif

    // 完全一致でループを抜け、max_lenに到達した場合
    return std::min(matched, max_len);
}
