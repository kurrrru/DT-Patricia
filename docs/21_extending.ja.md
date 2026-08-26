# DT-Patricia の拡張

DT-Patricia は 2 つのポリシーでパラメータ化されている。どの文字を区別するかを決める**アルファベット**と、各編集操作にいくらかかるかを決める**コストモデル**である。どちらも自分で書ける単なる構造体であり、ライブラリの内部に手を入れる必要はない。

このページは、それらの構造体が何を提供しなければならないかの規範的な仕様である。同梱のポリシーのどれを選ぶべきかについては [はじめに](10_getting_started.ja.md) を、同梱のポリシーの挙動については [API リファレンス](20_api_reference.ja.md) を参照。

## `AlphabetPolicy`

```cpp
#include <dt_patricia/policy/alphabet.hpp>

template <class Alphabet>
concept AlphabetPolicy = /* 後述 */;
```

アルファベットポリシーは、256 通りのバイト値すべてを `0 .. CODE_MAX` という詰まった範囲のコードへ写す対応を定める。ライブラリが文字を直接比較することはなく、比較するのはコードである。したがって同じコードに写る 2 文字は区別がつかない。`RyAlphabet` が `A` と `G` を畳み込むのも、各アルファベットがキャッチオールコードを持つのも、この仕組みによる。

木もアライナもポリシーでテンプレート化されているので、以上はすべてコンパイル時に解決され、テーブルはバイナリに焼き込まれる。

### 必要なメンバ

ポリシーは、ちょうど 4 つの public な静的メンバを持つクラスまたは構造体である。状態を持たず、実体化されることもない。

```cpp
static constexpr uint8_t CODE_TERM;
static constexpr uint8_t CODE_MAX;
static consteval std::array<uint8_t, 256> make_char_to_code() noexcept;
static consteval std::array<char, CODE_MAX + 1> make_code_to_char() noexcept;
```

- **`CODE_TERM`** — 終端記号のコード。**必ず 0 でなければならない。** 設定可能にせず固定しているのは、木のノードレイアウトが終端記号が先頭バケットを占めることに依存しているためである。
- **`CODE_MAX`** — ポリシーが生成しうる最大のコード。**0 より大きくなければならない。** 木はノードごとに `CODE_MAX + 1` 個の遷移スロットを確保するため、この値はメモリ使用量に直結する。例えば ASCII 値をそのままコードにするようなことは避け、コードの範囲は詰めて連続させること。
- **`make_char_to_code()`** — 符号化テーブルを返す。添字は `unsigned char` である。256 個の要素すべてを埋めなければならず、値はすべて `0 .. CODE_MAX` に収まらなければならない。大文字小文字・同義文字・未知の文字を畳み込むのはこのテーブルである。
- **`make_code_to_char()`** — 逆向きのテーブルを返す。各コードに対して、それを代表する 1 文字を与える。文字列の正規化とラベルの表示に使われる。

どちらのテーブル関数も `consteval` かつ `noexcept` であり、戻り値の型は上記のとおりでなければならない。コンセプトは `std::same_as` で戻り値型を検査するので、`auto` による型推論の結果配列の大きさが異なると、コンセプトを満たさない。

### 不変条件

コンセプトはメンバの存在確認で終わらない。コンパイル時に `detail::validate_alphabet_policy<Alphabet>()` も評価し、次を強制する。

1. `CODE_TERM == 0`。
2. `CODE_MAX > 0`。
3. `make_char_to_code()['\0'] == CODE_TERM`。
4. `make_char_to_code()` のすべての要素が `CODE_MAX` 以下。
5. `make_code_to_char()[CODE_TERM] == '\0'`。
6. **ラウンドトリップ整合性**: `0 .. CODE_MAX` のすべてのコード `c` について、`make_char_to_code()[make_code_to_char()[c]] == c`。

実際の間違いを捕まえるのは規則 6 である。これは、どのコードもそれ自身へ符号化し戻る代表文字を持たなければならない、と言っている。したがって特に、使われないコードがあってはならず、2 つのコードが代表文字を共有してもならない。一方のテーブルにコードを足してもう一方を忘れると、そのポリシーはコンセプトを満たさなくなり、それを用いた `PatriciaTree` の実体化はコンパイルに失敗する。

この検査は全体がコンパイル時に走るので、壊れたポリシーが実行時に誤った答えを返すことはありえない。単にビルドが通らないだけである。

### コンセプトが検査できない設計上の規則

- **未知の文字にキャッチオールコードを与えること。** テーブルは 256 通りのバイト値すべてを覆わなければならないので、名前を与えなかった文字も*どこか*へ行く必要がある。同梱のポリシーがそうしているように、それら専用のコードを 1 つ用意する。ただし帰結には注意すること。そうした文字は互いに一致するようになる。`ProteinAlphabet` のもとで `B` と `Z` が等しいのは、どちらもキャッチオールだからである。
- **`'\0'` 以外をコード 0 に写さないこと。** コード 0 は終端記号である。実在の文字がこれを共有すると木が壊れる。
- **文字ごとのコード定数を公開インターフェースに置かないこと。** 同梱のポリシーは `CODE_A` や `CODE_C` などを private に保っており、コンセプトも意図的にそれらを要求していない。private に保つことで、数値そのものは自由に変更できる状態が保たれる。
- **コードを詰めること。** `0 .. CODE_MAX` に穴があれば規則 6 でどのみち弾かれるが、より重要なのは、コードが 1 つ増えるごとにすべてのノードで遷移スロットが 1 つ増えるという点である。

### 実例

普通の英小文字テキスト向けのアルファベット。26 文字、大文字小文字を区別せず、それ以外はすべて畳み込む。

```cpp
#include <array>
#include <cstdint>

#include <dt_patricia/dt_patricia.hpp>

struct LowercaseAlphabet {
 public:
    static constexpr uint8_t CODE_TERM = 0;
    static constexpr uint8_t CODE_MAX = 27;  // 0 = 終端, 1..26 = 'a'..'z', 27 = その他

    static consteval std::array<uint8_t, 256> make_char_to_code() noexcept {
        std::array<uint8_t, 256> table{};
        table.fill(CODE_OTHER);

        table['\0'] = CODE_TERM;

        for (uint8_t i = 0; i < 26; ++i) {
            table[static_cast<unsigned char>('a' + i)] = static_cast<uint8_t>(i + 1);
            table[static_cast<unsigned char>('A' + i)] = static_cast<uint8_t>(i + 1);
        }

        return table;
    }

    static consteval std::array<char, CODE_MAX + 1> make_code_to_char() noexcept {
        std::array<char, CODE_MAX + 1> table{};
        table[CODE_TERM] = '\0';
        for (uint8_t i = 0; i < 26; ++i) {
            table[i + 1] = static_cast<char>('a' + i);
        }
        table[CODE_OTHER] = '?';  // その他すべての代表文字
        return table;
    }

 private:
    static constexpr uint8_t CODE_OTHER = 27;
};

static_assert(dt_patricia::AlphabetPolicy<LowercaseAlphabet>);
```

`static_assert` は必須ではないが、自分で書いたポリシーには必ず添えておくとよい。`PatriciaTree` の実体化の奥深くではなく、ポリシーの定義そのものの位置で違反が報告されるようになる。

使い方は同梱のポリシーと何ら変わらない。

```cpp
std::vector<std::string> dictionary = {"receive", "recieve", "retrieve"};
dt_patricia::PatriciaTree<LowercaseAlphabet> tree(dictionary);
dt_patricia::DTPatricia<LowercaseAlphabet, dt_patricia::UnitCost> aligner(tree);

auto hits = aligner.ed_within_k("recieve", 2);  // "receive" も見つかる
```

## コストポリシー

コストポリシーは意図的に、名前付きのコンセプトでは制約して**いない**。duck typing であり、どのメンバが必要かはポリシー自身が宣言するフラグによって変わる。

アルファベットポリシーと異なり、コストポリシーは実行時の状態を持つ。コスト値は通常の非静的メンバであり、アライナはそのコピーを保持する。

### 必要なメンバ

どのコストポリシーも、2 つのコンパイル時フラグを提供しなければならない。

```cpp
static constexpr bool is_linear;
static constexpr bool is_unit;
```

- **`is_linear`** — `L` 文字のギャップのコストが `gap * L` であるとき `true`。`false` は、ギャップを開くこと自体に追加コストがかかるアフィンモデルを選択する。このフラグは 2 つの `search_kernel` オーバーロード、すなわち別々の実装のどちらを使うかを決める。
- **`is_unit`** — `true` は専用の高速経路を選択する。**モデルがちょうど Levenshtein のとき、すなわち `is_linear == true` かつ `mismatch` と `gap` がともに 1 のときにのみ立てること。** 高速経路はコスト値を一切読まないため、それ以外の値と一緒に `is_unit = true` を宣言すると、指定したコストではなく単位コストの距離を黙って計算してしまう。ここに挙げた要件のうち、何もチェックしてくれないのはこれだけである。

フラグ以外に必要なメンバは、`is_linear` によって変わる。

| `is_linear` | 必要なメンバ | 長さ `L` のギャップのコスト |
| --- | --- | --- |
| `true` | `mismatch`, `gap` | `gap * L` |
| `false` | `mismatch`, `gap_open`, `gap_extend` | `gap_open + gap_extend * L` |

コストのメンバはいずれも `uint32_t` として読まれる。`UnitCost` のように `static constexpr` でもよいし、`LinearGapCost` や `AffineGapCost` のように構築時に設定される通常のメンバでもよい。実装は保持したコストオブジェクト越しに読むだけなので、どちらでも動く。

### 不変条件

これらを守るのは書き手の責任である。同梱のポリシーはコンストラクタで `std::invalid_argument` を投げて強制している。自分で書くポリシーも同様にするとよい。

- `mismatch >= 1`。
- `is_linear` のとき `gap >= 1`。枝刈りの段で `gap` による除算を行うため、0 は無意味であるにとどまらず致命的である。
- `!is_linear` のとき `gap_extend >= 1`。`gap_open` は 0 でもよく、その場合は線形モデルに退化する。
- `is_unit` ならば `mismatch == gap == 1`。

コストは、探索が保持する波面履歴の量も決める。保持されるレベル数は最大のコスト値とともに増える。非常に大きなコスト値も正しく動くが、それに比例してメモリを消費する。望む比を表現できる最小の整数を選ぶこと。`LinearGapCost(100, 300)` ではなく `LinearGapCost(1, 3)` とする。

### 実例

コスト比をコンパイル時に固定した線形モデル。すべて `static constexpr` なので、この型はデフォルト構築可能であり、コストオブジェクトを渡さずに使える。

```cpp
struct GapHeavyCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit = false;  // unit ではない: コストが 1/1 ではないため
    static constexpr uint32_t mismatch = 1;
    static constexpr uint32_t gap = 3;
};

dt_patricia::DTPatricia<dt_patricia::DnaAlphabet, GapHeavyCost> aligner(tree);
```

実行時に設定するものは `LinearGapCost` と同じ形になる。

```cpp
struct MyCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit = false;
    uint32_t mismatch;
    uint32_t gap;

    MyCost(uint32_t m, uint32_t g) : mismatch(m), gap(g) {
        if (m < 1 || g < 1) {
            throw std::invalid_argument("MyCost: mismatch and gap must be >= 1");
        }
    }
};

dt_patricia::DTPatricia<dt_patricia::DnaAlphabet, MyCost> aligner(tree, MyCost(2, 5));
```

この型はデフォルトコンストラクタを持たないので、コスト引数は必須になる。アライナ側の既定引数は、省略した場合にのみ実体化されるためである。

## 新しいポリシーの検証

このライブラリ自身のテストは、結果を素朴な動的計画法の実装と突き合わせて検証している。自分で書いたポリシーに確信を持つには、同じやり方をするのがよい。対象のアルファベット上でランダムな辞書とクエリを生成し、同じコスト値を使った素朴な `O(nm)` の DP で距離を計算し、両者が厳密に一致することを確認する。既存の枠組みは `tests/` にある。

アルファベットポリシーに関しては、`canonicalize<Alphabet>` が最も手軽な健全性チェックになる。ライブラリが実際に何を比較するかがそのまま見えるので、`canonicalize<LowercaseAlphabet>("Hello!")` が `"hello?"` を返すことを見れば、畳み込みの規則を一行で確認できる。
