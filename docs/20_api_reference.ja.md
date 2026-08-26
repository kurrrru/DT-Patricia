# API リファレンス

DT-Patricia の公開インターフェースの正確な契約。何をするものかではなく、どうやって目的を達成するかを知りたい場合は [はじめに](10_getting_started.ja.md) から読むこと。

以下はすべて名前空間 `dt_patricia` に属する。

## 目次

- [ヘッダ](#ヘッダ)
- [`AlignmentResult`](#alignmentresult)
- [`PatriciaTree<Alphabet>`](#patriciatreealphabet)
- [`DTPatricia<Alphabet, CostType>`](#dtpatriciaalphabet-costtype)
- [クエリ結果のセマンティクス](#クエリ結果のセマンティクス)
- [コストポリシー](#コストポリシー)
- [アルファベットポリシー](#アルファベットポリシー)
- [自由関数](#自由関数)
- [スレッド安全性](#スレッド安全性)
- [例外](#例外)

## ヘッダ

| ヘッダ | 内容 |
| --- | --- |
| `<dt_patricia/dt_patricia.hpp>` | エントリポイント。以下すべてをインクルードする。通常はこれを使う |
| `<dt_patricia/patricia_tree.hpp>` | `PatriciaTree` |
| `<dt_patricia/aligner.hpp>` | `DTPatricia` |
| `<dt_patricia/alignment_result.hpp>` | `AlignmentResult` |
| `<dt_patricia/policy/alphabet.hpp>` | `AlphabetPolicy`、同梱のアルファベット、`canonicalize` |
| `<dt_patricia/policy/cost.hpp>` | `UnitCost`、`LinearGapCost`、`AffineGapCost` |

`dt_patricia/internal/` と `dt_patricia/debug/` は公開インターフェースではない。その内容は予告なく変わりうるため、これらをインクルードしたり、`dt_patricia::internal` の名前を参照したりすることは推奨しない。

## `AlignmentResult`

```cpp
struct AlignmentResult {
    uint32_t string_id;
    uint32_t score;
};
```

- **`string_id`** — 一致した文字列の、`PatriciaTree` のコンストラクタに渡した `std::vector<std::string>` における（入力順の）インデックス。内部のノード ID ではない。
- **`score`** — クエリとその文字列との、アライナのコストモデルのもとでの厳密な最適アラインメントコスト。`UnitCost` の場合は Levenshtein 距離になる。

集成体であり、メンバ関数を持たない。

## `PatriciaTree<Alphabet>`

```cpp
template <AlphabetPolicy Alphabet = DnaAlphabet>
class PatriciaTree;
```

辞書を保持する、変更不可能なパス圧縮トライ。まずこれを構築し、その上に 1 つ以上のアライナを作る。

`Alphabet` は [`AlphabetPolicy`](21_extending.ja.md#alphabetpolicy) コンセプトを満たさなければならない。

### 構築と特殊メンバ

```cpp
explicit PatriciaTree(const std::vector<std::string> &input_data);

PatriciaTree() = delete;
PatriciaTree(const PatriciaTree &) = delete;             // 辞書は大きくなりうるため
PatriciaTree &operator=(const PatriciaTree &) = delete;
PatriciaTree(PatriciaTree &&) noexcept = default;
PatriciaTree &operator=(PatriciaTree &&) noexcept = default;
~PatriciaTree() = default;
```

コンストラクタが木全体を構築する。別途ビルド手順はなく、構築後に挿入・削除する手段もない。`input_data` は保持されず（木は必要なものを複製する）、コンストラクタから戻った時点で破棄してよい。ただし string ID は `input_data` 内の位置を指すので、文字列そのものが必要な呼び出し側は自分で保持しておく必要がある。

`input_data` は空でもよく、空文字列を含んでもよく、重複を含んでもよい。重複は別々の ID を保ち、別々の結果として報告される。

どの文字列も `'\0'` を含んではならない。

木はムーブ可能だがコピー不可能である。

### メンバ型と定数

```cpp
using alphabet_type = Alphabet;

static constexpr uint8_t     CODE_TERM = Alphabet::CODE_TERM;   // 常に 0
static constexpr uint8_t     CODE_MAX  = Alphabet::CODE_MAX;
static constexpr std::size_t BUCKET_SIZE = CODE_MAX + 1;
static constexpr std::size_t SIMD_PADDING_SIZE = 256;

inline static constexpr std::array<uint8_t, 256> CHAR_TO_CODE = Alphabet::make_char_to_code();
```

`SIMD_PADDING_SIZE` は、ベクトル化されたロードが末尾を越えて読まないように、探索カーネルが内部に持つクエリのコピーへ付加する余白の大きさである。これは呼び出し側への要求ではなく定数として公開されている。渡すクエリ文字列に自前のパディングは不要である。

### 辞書レベルの問い合わせ

```cpp
[[nodiscard]] bool     empty() const noexcept;
[[nodiscard]] uint32_t string_count() const noexcept;
[[nodiscard]] uint32_t node_count() const noexcept;
```

- **`empty()`** — 木がノードを 1 つも持たない場合、すなわち空のベクタから構築された場合に `true`。空の木へのクエリは定義された動作であり、結果を返さない。
- **`string_count()`** — 辞書中の文字列数。重複も別々に数える。`input_data.size()` に等しい。
- **`node_count()`** — 予約された ID 0 を含む、木のノードスロット数。辞書の大きさを表す意味のある指標ではなく、レイアウトの実装詳細である。

### ノードレベルの走査

これらは木の構造そのものを露出する。辞書の検査・エクスポート・独自の走査の実装に有用なため公開しているが、ライブラリの通常の利用者が必要とすることはない。木を利用する意図された手段は `DTPatricia` の探索である。

ノード ID `0` は予約された無効値であり、根は ID `1` である。以下のアクセサはいずれも有効なノード ID を前提とする。注記のある場合を除き、`0` や範囲外の ID を渡すのは未定義動作である。

```cpp
[[nodiscard]] uint32_t root_id() const noexcept;                              // 常に 1
[[nodiscard]] uint32_t transition(uint32_t node_id, char ch) const noexcept;
[[nodiscard]] uint32_t transition(uint32_t node_id, uint8_t code) const noexcept;
[[nodiscard]] uint32_t get_parent(uint32_t node_id) const noexcept;
[[nodiscard]] bool     is_leaf(uint32_t node_id) const noexcept;
[[nodiscard]] bool     is_terminal(uint32_t node_id) const noexcept;
[[nodiscard]] std::string_view get_label(uint32_t node_id) const noexcept;
[[nodiscard]] uint32_t get_label_length(uint32_t node_id) const noexcept;
[[nodiscard]] std::span<const uint32_t> get_string_id(uint32_t node_id) const noexcept;
```

- **`transition(node_id, ch)`** — `node_id` から出る辺のうち先頭文字が `ch` であるものをたどり、到達するノードの ID を返す。そのような辺がなければ `0` を返す。`char` を取るオーバーロードは `ch` をまず `CHAR_TO_CODE` で符号化するので、アルファベットの畳み込み規則（大文字小文字、`U`/`T`、キャッチオール）が適用される。`uint8_t` を取るオーバーロードは、符号化済みのコードを受け取る。
- **`get_parent(node_id)`** — 親のノード ID。根および範囲外の ID に対しては `0` を返す。
- **`is_leaf(node_id)`** — そのノードが出る辺を持たない場合に `true`。
- **`is_terminal(node_id)`** — 辞書中の文字列がこのノードで終わる場合、すなわち `CODE_TERM` の出辺を持つ場合に `true`。ただし string ID は終端ノード自身ではなく、その辺の先の葉に紐づいている点に注意する（`get_string_id(transition(node_id, CODE_TERM))`）。
- **`get_label(node_id)`** — `node_id` に**入ってくる**辺のラベル。返されるビューは木が保持する記憶域を指し、木が生存していてムーブ元になっていない限り有効である。
- **`get_string_id(node_id)`** — この葉で終わる辞書文字列の ID 群。`node_id` が葉でない場合は空の span を返す。辞書に重複があった場合は複数の ID が返る。

### 事前計算された部分木の統計

```cpp
[[nodiscard]] const std::vector<uint32_t> &get_subtree_counts() const noexcept;
[[nodiscard]] const std::vector<uint32_t> &get_subtree_max_lengths() const noexcept;
[[nodiscard]] const std::vector<uint32_t> &get_subtree_min_lengths() const noexcept;
[[nodiscard]] const std::vector<uint32_t> &get_parent_path_lengths() const noexcept;
```

いずれもノード ID を添字とし、構築時に一度だけ計算される。探索はこれらを枝刈りに用いる。公開している理由は上記の走査アクセサと同じである。

- **`get_subtree_counts()`** — 各ノードを根とする部分木に含まれる辞書文字列の数。
- **`get_subtree_max_lengths()` / `get_subtree_min_lengths()`** — 各ノードから、その部分木の葉までの距離の最大値と最小値。そのノード自身の入辺ラベルを含む。
- **`get_parent_path_lengths()`** — 根から各ノードの**親**までのパス長。そのノード自身の入辺ラベルは含まない。

返される参照は、木が生存していてムーブ元になっていない限り有効である。

## `DTPatricia<Alphabet, CostType>`

```cpp
template <AlphabetPolicy Alphabet = DnaAlphabet, typename CostType = UnitCost>
class DTPatricia;
```

探索エンジン。`PatriciaTree` の上で diagonal transition アルゴリズムを実行する。

`Alphabet` は木を実体化したものと**同じ**型でなければならない。食い違いはコンパイルエラーになる。`CostType` は[コストポリシーの要件](21_extending.ja.md#コストポリシー)を満たさなければならない。

### 構築と特殊メンバ

```cpp
DTPatricia(const tree_type &patricia_tree, CostType cost = CostType());

DTPatricia() = delete;
DTPatricia(const DTPatricia &) = delete;
DTPatricia &operator=(const DTPatricia &) = delete;
DTPatricia(DTPatricia &&) noexcept = delete;
DTPatricia &operator=(DTPatricia &&) noexcept = delete;
~DTPatricia() = default;
```

アライナは木への**参照**と、コストオブジェクトの**コピー**を保持する。木はアライナより長く生存していなければならない。構築は軽量であり、アライナごとのインデックスは作られない。したがって、1 つの木の上に異なるコストモデルのアライナを複数作ることが可能である。

その参照を保持しているがゆえに、コピーもムーブもできない。

`cost` の既定引数は、引数を省略した場合にのみ実体化される。したがって、デフォルトコンストラクタを持たない `CostType`（`LinearGapCost` と `AffineGapCost` はいずれもそうである）は、単にコストを渡すことを要求するだけである。

### メンバ型

```cpp
using alphabet_type = Alphabet;
using tree_type     = PatriciaTree<Alphabet>;
```

### `get_patricia_tree`

```cpp
[[nodiscard]] const tree_type &get_patricia_tree() const noexcept;
```

アライナが構築対象とした木。

### `ed_to_all`

```cpp
std::vector<AlignmentResult> ed_to_all(const std::string &query) const;
```

`query` から辞書中の**すべての**文字列までの距離を返す。結果はちょうど `string_count()` 件であり、木が空なら 0 件である。

枝刈りに使える閾値がないため、3 つのクエリの中で最も高コストである。

### `ed_within_k`

```cpp
std::vector<AlignmentResult> ed_within_k(const std::string &query, int k) const;
```

`query` からの距離が **`k` 以下**（境界値を含む）である辞書文字列をすべて返す。`k` は編集操作の回数ではなくアライナのコストモデルにおけるコストなので、`LinearGapCost(1, 3)` のもとではギャップ 1 つですでにコスト 3 である。

`k` は枝刈りの上界も兼ねており、その枝刈りは `k` に関して厳密である。閾値以内の文字列が捨てられることはない。したがって `k` が小さいほど探索は速くなる。

停止条件が評価されるのは常にある距離レベルを完了した後なので、`k == 0` は完全一致を返し、`k < 0` は何も返さないのではなく `k == 0` と同じ挙動になる。

### `ed_kth_smallest`

```cpp
std::vector<AlignmentResult> ed_kth_smallest(const std::string &query, size_t k) const;
```

最も近い `k` 件の辞書文字列を返す。`k` は `string_count()` に切り詰められる。上と同じ理由により、`k == 0` は何も返さず終わるのではなく完全一致を返す。

**結果が `k` 件を超えることがある。** 探索は、十分な件数が集まったかを判定できるようになる前に各距離レベルを完了させる。そのため、打ち切り距離で複数の文字列が同点になった場合、それらはすべて返される。これは、等しく良い候補の間で恣意的な選択をしないためである。ちょうど `k` 件が必要なら、呼び出し側で結果を切り詰めること。

打ち切り距離は事前にはわからないため、枝刈りの上界は使われない。

### `search_kernel`

```cpp
template <typename StopPredicate>
std::vector<AlignmentResult> search_kernel(const std::string &query,
                                           StopPredicate stop_predicate,
                                           int upper_bound = -1) const;
```

他の 3 つのクエリが薄くラップしているエンジン本体。それらでは表現できない停止規則が必要なときに使う。

- **`stop_predicate`** — `bool(int current_score, const std::vector<AlignmentResult> &results)` というシグネチャで呼び出せるもの。**距離レベルごとに 1 回**、その距離のすべての一致が `results` に追加された後で呼ばれる。`true` を返すと探索を打ち切り、それまでに集まった結果がそのまま返される。`current_score` は処理し終えた距離レベルであり、呼び出しをまたいで単調に増加する。レベルの合間にしか参照されないため、これを通じて表現される予算——実時間の締め切りを含む——の粒度は「距離レベル 1 つ分の仕事」になる。
- **`upper_bound`** — `ed_within_k` と同様に、コスト `upper_bound` での枝刈りを有効にする。既定値の `-1` は枝刈りを無効にする。枝刈りが上界以内でありうる候補を捨てることはない。ただし、枝刈り自体は探索を停止させ**ない**点に注意する。上界を尊重する `stop_predicate` も併せて渡さないと、探索は見つけるものがなくなった後も上界を越えて進み続ける。

`CostType::is_linear` で制約された 2 つのオーバーロードが存在するが、適切なほうが自動的に選ばれ、呼び出し側からは区別できない。

タイムアウトを表す述語と、top-k と閾値を組み合わせた述語の実例は [`examples/basic_example.cpp`](../examples/basic_example.cpp) の 6 節と 7 節にある。

## クエリ結果のセマンティクス

以下は 4 つのクエリ関数すべてに当てはまる。

- **順序**: 結果は `score` の非減少順に追加される。同じ score を持つエントリ同士の相対順序は未規定であり、依存してはならない。
- **一意性**: 各 `string_id` は高々 1 回しか現れない。辞書中の重複した文字列は異なる ID を持つので、それぞれが 1 回ずつ現れる。
- **厳密性**: 報告される `score` はすべて最適コストである。近似もヒューリスティックも用いない。
- **空の木**: 空のベクタを返す。
- **空のクエリ**: 定義された動作である。各辞書文字列までの距離は、その文字列全体を削除するコストになる。
- **クエリ文字列**: クエリは内部でコピーされ正規化される。渡した文字列が変更されることはない。

## コストポリシー

`<dt_patricia/policy/cost.hpp>` で定義される。コストはすべて `uint32_t` である。

### `UnitCost`

```cpp
struct UnitCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit   = true;
    static constexpr uint32_t mismatch = 1;
    static constexpr uint32_t gap      = 1;
};
```

素の Levenshtein 距離。状態を持たずデフォルト構築可能なので、`DTPatricia<Alphabet, UnitCost>` はコスト引数を必要としない。`is_unit` が専用の実行経路を選択するため、同じコストを `LinearGapCost` で表現するより明確に速い。

### `LinearGapCost`

```cpp
struct LinearGapCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit   = false;
    uint32_t mismatch;
    uint32_t gap;

    LinearGapCost(uint32_t m, uint32_t g);
};
```

`L` 文字のギャップのコストは `gap * L` である。`m < 1` または `g < 1` の場合 `std::invalid_argument` を投げる。デフォルト構築はできない。

### `AffineGapCost`

```cpp
struct AffineGapCost {
    static constexpr bool is_linear = false;
    static constexpr bool is_unit   = false;
    uint32_t mismatch;
    uint32_t gap_open;
    uint32_t gap_extend;

    AffineGapCost(uint32_t m, uint32_t g_open, uint32_t g_extend);
};
```

`L` 文字のギャップのコストは `gap_open + gap_extend * L` である。`m < 1` または `g_extend < 1` の場合 `std::invalid_argument` を投げる。`gap_open` は 0 でもよい。デフォルト構築はできない。

## アルファベットポリシー

`<dt_patricia/policy/alphabet.hpp>` で定義される。いずれも 256 通りのバイト値すべてを小さなコードの集合へ写す。コード `0` は常に終端記号である。同じコードを持つ文字は探索から見て区別がつかない。また各アルファベットはキャッチオールコードを 1 つ持ち、名前を与えられていない文字はすべてそこへ落ちる。したがって、それらの文字は互いに一致する。

大文字と小文字はいずれのアルファベットでも区別されない。

| ポリシー | `CODE_MAX` | 名前を与えられた文字 | キャッチオールの代表文字 |
| --- | --- | --- | --- |
| `DnaAlphabet` | 5 | `A` `C` `G` `T`。`U` は `T` に畳まれる | `N` |
| `RyAlphabet` | 3 | `R`（`A`, `G`, `R` から）、`Y`（`C`, `T`, `U`, `Y` から） | `N` |
| `ProteinAlphabet` | 21 | `A C D E F G H I K L M N P Q R S T V W Y` | `X` |

`ProteinAlphabet` のもとでは、曖昧記号 `B`、`J`、`O`、`U`、`X`、`Z` はすべてキャッチオールに落ちるため互いに等しい。

これらが満たす `AlphabetPolicy` コンセプトについては [DT-Patricia の拡張](21_extending.ja.md#alphabetpolicy) を参照。

## 自由関数

```cpp
template <AlphabetPolicy Alphabet>
std::string canonicalize(std::string_view s);

template <AlphabetPolicy Alphabet>
void canonicalize_inplace(char *data, std::size_t n) noexcept;
```

各文字を、`Alphabet` のもとでのそのコードの代表文字に書き換える。探索がクエリに対して内部で行っているのがこれであり、ライブラリが実際に何を比較しているかを呼び出し側が確認できるように公開されている。`canonicalize<DnaAlphabet>("acgu")` は `"ACGT"` であり、`RyAlphabet` のもとでは `"RYRY"` になる。

2 つの文字列の距離が 0 であることと、正規形が等しいことは同値である。

`canonicalize_inplace` は `data` から `n` バイトを上書きする。変換表は静的記憶域期間を持つコンパイル時定数なので、どちらの関数も複数スレッドから安全に呼べる。

## スレッド安全性

`PatriciaTree` は構築後に変更されない。また `DTPatricia` のクエリ関数はいずれも `const` であり、作業状態をすべてローカル変数に持つ。したがって、同じ木、さらには同じアライナに対して、外部同期なしに任意個のスレッドから同時にクエリを投げてよい。

内部での並列化は行っていない。1 つのクエリは 1 つのスレッドで実行される。

## 例外

ライブラリが例外を投げるのは、正でないコストに対して `std::invalid_argument` を送出するコストポリシーのコンストラクタだけである。標準ライブラリからの確保失敗は通常どおり伝播する。クエリ関数自体は例外を投げない。
