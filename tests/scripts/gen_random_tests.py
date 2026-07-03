#!/usr/bin/env python3
# ============================================================
# gen_random_tests.py
#   シード固定の擬似乱数で brute_test 用のランダムテストケースを量産する。
#
#   生成物は tests/testcase_random/*.txt。test.cpp の想定するフォーマット
#   （ALPHABET / COST / TARGETS / QUERIES / TESTS）で出力する。
#   これらを brute_test の引数に testcase_random ディレクトリを渡すことで、
#   DT-Patricia と愚直 DP の一致をランダムケースでも検証できる。
#
# [使い方]
#   # デフォルト（シード 42, 既定件数）で生成
#   python3 tests/scripts/gen_random_tests.py
#
#   # シードや件数、出力先を指定
#   python3 tests/scripts/gen_random_tests.py --seed 42 --count 80
#
#   # 生成後にランダムテストを実行（tests/ をカレントにして）
#   ./build/tests/tests/brute_test testcase_random
#
# [方針]
#   - test.cpp の比較相手は O(n*m) の愚直 DP。配列長が長いと計算が終わらないため、
#     配列長は最大でも 100 程度に抑える。長いケースはターゲット数・クエリ数を減らし、
#     短いケースはターゲットやクエリを多くして「複数ターゲット・大量クエリ」を量産する。
#   - すべての乱数は単一の random.Random(seed) から順番に消費するので、
#     同じシード・同じ件数なら生成物は完全に再現する。
# ============================================================

import argparse
import os
import random

# ------------------------------------------------------------
# アルファベットごとの使用可能文字（test.cpp の allowed_chars と一致させる）
#   DnaAlphabet     : A/C/G/T
#   RyAlphabet      : A/C/G/T/U/R/Y（A,G,R と C,T,U,Y がそれぞれ同一コードに潰れる）
#   ProteinAlphabet : 20 アミノ酸
# ------------------------------------------------------------
ALPHABETS = {
    "DnaAlphabet": "ACGT",
    "RyAlphabet": "ACGTURY",
    "ProteinAlphabet": "ACDEFGHIKLMNPQRSTVWY",
}
# アルファベットの出現重み（DNA を中心に、Ry / Protein も混ぜる）
ALPHABET_WEIGHTS = {
    "DnaAlphabet": 6,
    "RyAlphabet": 2,
    "ProteinAlphabet": 2,
}

# ------------------------------------------------------------
# 配列長の上限（愚直 DP を現実的な時間で終わらせるための制約）
# ------------------------------------------------------------
MAX_LEN = 100


# ------------------------------------------------------------
# ケースの「形」。長さと個数のトレードオフでカテゴリを分ける。
#   nt : ターゲット数の範囲
#   nq : クエリ数の範囲
#   ln : 各配列長の範囲（上限は MAX_LEN で最終的にクランプ）
# ------------------------------------------------------------
SHAPES = [
    # 短い配列を大量のターゲットで（複数ターゲットの検証）
    {"name": "many_targets_short", "nt": (20, 60), "nq": (3, 12), "ln": (3, 20)},
    # 短い配列を大量のクエリで（大量クエリの検証）
    {"name": "many_queries_short", "nt": (3, 12), "nq": (20, 50), "ln": (3, 20)},
    # 中規模を程よい個数で
    {"name": "medium_mixed", "nt": (5, 20), "nq": (3, 10), "ln": (10, 50)},
    # 長い配列は個数を絞る（DP が重いため）
    {"name": "long_few", "nt": (1, 6), "nq": (1, 4), "ln": (50, MAX_LEN)},
    # 長さ 1〜5 の極端に短いエッジケース
    {"name": "tiny_edge", "nt": (2, 10), "nq": (2, 8), "ln": (1, 5)},
]


# ------------------------------------------------------------
# コスト設定をランダムに選ぶ。
#   test.cpp / cost.hpp の制約:
#     LinearGapCost : mismatch>=1, gap>=1
#     AffineGapCost : mismatch>=1, gap_extend>=1（gap_open は 0 も許容）
#   愚直 DP と DT-Patricia は同じコスト式を使うので、この範囲なら比較は妥当。
#   極端な非対称コストで偽陽性を出さないよう、常識的な範囲に収める。
# ------------------------------------------------------------
def pick_cost(rng):
    kind = rng.choices(["UnitCost", "LinearGapCost", "AffineGapCost"],
                       weights=[3, 3, 3])[0]
    if kind == "UnitCost":
        return "COST UnitCost"
    if kind == "LinearGapCost":
        gap = rng.randint(1, 3)
        # mismatch <= 2*gap にしておくと置換・ギャップ2連の両方が現れやすい
        mismatch = rng.randint(1, 2 * gap)
        return f"COST LinearGapCost {mismatch} {gap}"
    # AffineGapCost: gap_open + gap_extend * L
    extend = rng.randint(1, 3)
    gap_open = rng.randint(0, 4)
    mismatch = rng.randint(1, gap_open + 2 * extend)
    return f"COST AffineGapCost {mismatch} {gap_open} {extend}"


# ------------------------------------------------------------
# 配列生成ヘルパ
# ------------------------------------------------------------
def rand_len(rng, lo, hi):
    """[lo, hi] の長さを MAX_LEN 以内・1 以上でクランプして返す。"""
    return max(1, min(MAX_LEN, rng.randint(lo, hi)))


def rand_seq(rng, chars, length):
    return "".join(rng.choice(chars) for _ in range(length))


# ------------------------------------------------------------
# TESTS セクション（実行する操作）を生成する。
#   ed_to_all は必ず含める。ed_within_k / ed_kth_smallest は k を散らして加える。
#   op 数はランタイム抑制のため 5 個までに制限する。
# ------------------------------------------------------------
def gen_ops(rng, num_targets, max_len):
    ops = []
    # within_k: 0（完全一致のみ）から大きめまで散らす
    within_ks = {0, 1, rng.randint(1, max(1, max_len // 4)), rng.randint(1, max(1, max_len))}
    for k in sorted(within_ks):
        if rng.random() < 0.6:
            ops.append(f"ed_within_k {k}")
    # kth_smallest: 1、中間、件数ちょうど、件数超過（クランプ確認）
    kth_ks = {1, max(1, num_targets // 2), num_targets, num_targets + rng.randint(1, 2)}
    for k in sorted(kth_ks):
        if rng.random() < 0.5:
            ops.append(f"ed_kth_smallest {k}")

    rng.shuffle(ops)
    ops = ops[:4]  # ed_to_all を足して最大 5 個
    ops.append("ed_to_all")
    return ops


# ------------------------------------------------------------
# 1 ファイル分のテキストを組み立てる
# ------------------------------------------------------------
def build_case(rng, index):
    shape = rng.choice(SHAPES)
    alphabet = rng.choices(list(ALPHABET_WEIGHTS), weights=list(ALPHABET_WEIGHTS.values()))[0]
    chars = ALPHABETS[alphabet]
    cost_line = pick_cost(rng)

    # カテゴリの個数・長さ範囲に従い、ターゲット・クエリをそれぞれ独立に生成する
    nt = rng.randint(*shape["nt"])
    nq = rng.randint(*shape["nq"])
    lo, hi = shape["ln"]
    targets = [rand_seq(rng, chars, rand_len(rng, lo, hi)) for _ in range(nt)]
    queries = [rand_seq(rng, chars, rand_len(rng, lo, hi)) for _ in range(nq)]
    max_len = max(len(s) for s in targets + queries)
    ops = gen_ops(rng, len(targets), max_len)

    lines = []
    lines.append(f"# 自動生成ランダムテスト #{index:03d} (shape={shape['name']})")
    lines.append(f"# generated by tests/scripts/gen_random_tests.py")
    lines.append(f"ALPHABET {alphabet}")
    lines.append(cost_line)
    lines.append("")
    lines.append("TARGETS")
    lines.extend(targets)
    lines.append("")
    lines.append("QUERIES")
    lines.extend(queries)
    lines.append("")
    lines.append("TESTS")
    lines.extend(ops)
    lines.append("")
    return "\n".join(lines)


# ------------------------------------------------------------
# main
# ------------------------------------------------------------
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_out = os.path.normpath(os.path.join(script_dir, "..", "testcase_random"))

    parser = argparse.ArgumentParser(description="brute_test 用ランダムテストケース生成")
    parser.add_argument("--seed", type=int, default=42,
                        help="乱数シード（固定して再現性を担保する）")
    parser.add_argument("--count", type=int, default=60,
                        help="生成するテストファイル数")
    parser.add_argument("--out-dir", default=default_out,
                        help="出力先ディレクトリ（既定: tests/testcase_random）")
    parser.add_argument("--clean", action="store_true", default=True,
                        help="出力先の既存 rand_*.txt を先に削除する（既定: 有効）")
    parser.add_argument("--no-clean", dest="clean", action="store_false")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    if args.clean:
        for name in os.listdir(args.out_dir):
            if name.startswith("rand_") and name.endswith(".txt"):
                os.remove(os.path.join(args.out_dir, name))

    rng = random.Random(args.seed)
    for i in range(args.count):
        text = build_case(rng, i)
        path = os.path.join(args.out_dir, f"rand_{i:03d}.txt")
        with open(path, "w") as f:
            f.write(text)

    print(f"generated {args.count} random test files -> {args.out_dir} (seed={args.seed})")


if __name__ == "__main__":
    main()
