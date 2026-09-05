# はじめに

このページは、空のプロジェクトから DT-Patricia のクエリが動くところまでを案内する。上から順に読むこと。掲載したコードはすべて、そのままコンパイルして実行できる。

ここで触れた要素の正確なシグネチャと挙動については [API リファレンス](20_api_reference.ja.md) を参照。独自のアルファベットやコストモデルを組み込みたい場合は [DT-Patricia の拡張](21_extending.ja.md) を参照。

## 1. 要件

- C++20 に対応したコンパイラ。CI では Ubuntu 上の GCC と Clang を検証している。
- CMake 3.21 以降（`CMakePresets.json` を使う場合は 3.25 以降）。

DT-Patricia はヘッダオンリーで、外部ライブラリへの依存はない。リンクすべきバイナリは生成されないので、ヘッダをインクルードパスに追加するだけでよい。

## 2. プロジェクトへの導入

依存管理の方法に合わせて、以下のいずれかを選ぶ。

### 方法 A: ソースを同梱する（`add_subdirectory`）

リポジトリをプロジェクト内の適当な場所（例えば `external/DT-Patricia`）に置き、`CMakeLists.txt` から追加する。

```cmake
add_subdirectory(external/DT-Patricia)
target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

`dt_patricia::dt_patricia` ターゲットにリンクすると、インクルードパスと C++20 の要求の両方が伝播する。`CMAKE_CXX_STANDARD` を自分で設定する必要はない。

この方法で取り込んだ場合、DT-Patricia はトップレベルのプロジェクトではなくなるため、examples と tests はビルドされない。

### 方法 B: configure 時に取得する（`FetchContent`）

```cmake
include(FetchContent)
FetchContent_Declare(dt_patricia
    GIT_REPOSITORY https://github.com/kurrrru/DT-Patricia.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(dt_patricia)

target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

### 方法 C: インストールして `find_package` する

ヘッダと CMake パッケージファイルを一度インストールする。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/prefix
cmake --install build
```

以降、任意のプロジェクトから次のように使える。

```cmake
find_package(dt_patricia REQUIRED)
target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

インストール先が CMake の既定の探索対象でない場合は、利用側プロジェクトの configure 時に指定する。

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/prefix
```

インストールされるパッケージのバージョン互換性は `SameMajorVersion` なので、`find_package(dt_patricia 0.1 REQUIRED)` は任意の `0.x` のインストールを受け入れる。

### 方法 D: CMake を使わない

`include/` をコンパイラのインクルードパスに追加し、C++20 でコンパイルする。

```bash
g++ -std=c++20 -O2 -Iexternal/DT-Patricia/include my_app.cpp -o my_app
```

## 3. 最初のクエリ

エントリポイントとなるヘッダは一つだけである。これをインクルードすれば、このページで扱う内容はすべて使える。

```cpp
#include <dt_patricia/dt_patricia.hpp>
```

このライブラリの使い方は、常に次の 3 ステップである。

1. 辞書（`std::vector<std::string>`）から `PatriciaTree` を構築する。
2. その木の上に `DTPatricia` アライナを作る。
3. アライナに対してクエリを実行する。

```cpp
#include <iostream>
#include <string>
#include <vector>

#include <dt_patricia/dt_patricia.hpp>

int main() {
    using namespace dt_patricia;

    std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG", "TGCA"};

    // 1) 辞書のインデックスを構築する
    PatriciaTree<DnaAlphabet> tree(targets);

    // 2) その上にアライナを作る
    DTPatricia<DnaAlphabet, UnitCost> aligner(tree);

    // 3) クエリ: "ACGT" から編集距離 1 以内の辞書エントリをすべて取得する
    for (const AlignmentResult &r : aligner.ed_within_k("ACGT", 1)) {
        std::cout << r.score << "  " << targets[r.string_id] << "\n";
    }
}
```

期待される出力（距離 1 の 3 件は、その 3 件の間では任意の順序で並びうる）。

```
0  ACGT
1  AAGT
1  ACG
1  ACGA
```

結果型については、この時点で次の 2 点を押さえておくとよい。

- `string_id` は、**木に渡したベクタへのインデックス**であり、渡したときの順序に対応する。ライブラリが文字列そのものを返すことはないので、引くのは呼び出し側の責任である。
- `score` は、選んだコストモデルのもとでの**厳密な**編集距離である。推定値ではなく、閾値以内の候補を取りこぼすこともない。

結果は `score` の非減少順に返り、各辞書エントリは高々 1 回しか現れない。

### ビルドと実行

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/my_app
```

`Release`（少なくとも最適化を有効にした構成）でビルドすること。探索カーネルはコンパイラがベクトル化できるように書かれており、最適化なしのビルドは著しく遅くなる。

## 4. アルファベットを選ぶ

アルファベットはコンパイル時のポリシーであり、**どの文字を区別するか**を決める。同じコードに写る 2 文字は等しいものとして扱われるため、アルファベットを変えることは距離の意味そのものを変えることになる。

ライブラリには 3 つのアルファベットが同梱されている。

| アルファベット | 区別するもの | 使いどころ |
| --- | --- | --- |
| `DnaAlphabet` | `A` `C` `G` `T` と、その他をまとめる 1 コード | DNA/RNA 配列。`U` は `T` に畳まれる |
| `RyAlphabet` | プリン (`A`,`G`,`R`) とピリミジン (`C`,`T`,`U`,`Y`)、およびその他 | プリン/ピリミジンのパターンだけを見たい場合 |
| `ProteinAlphabet` | 標準的な 20 種のアミノ酸と、その他をまとめる 1 コード | タンパク質配列 |

核酸なら `DnaAlphabet`、ペプチドなら `ProteinAlphabet` から始めるとよい。`RyAlphabet` は、2 クラスへの縮約が本当に目的である場合にのみ選ぶ。このアルファベットのもとでは `GCAT` も `ACGT` もともに `RYRY` であり、距離は 0 になる。

いずれも当てはまらない場合、ポリシークラスを 1 つ書くだけで任意のアルファベットに対応できる。[DT-Patricia の拡張](21_extending.ja.md) を参照。

なお、アルファベットは型の一部である。`DTPatricia<Alphabet, Cost>` は、**同じ** `Alphabet` で実体化された `PatriciaTree<Alphabet>` の上にしか構築できない。食い違いは黙って誤った答えを返すのではなく、コンパイルエラーになる。

## 5. コストモデルを選ぶ

コストモデルは、**各編集操作にいくらかかるか**を決める。3 つが提供されている。

`UnitCost` — mismatch と gap がともにコスト 1、すなわち素の Levenshtein 距離。既定であり、かつ最も高速でもある。他の 2 つでは使えない専用の実行経路を通るためである。特に理由がなければこれを使う。

```cpp
DTPatricia<DnaAlphabet, UnitCost> aligner(tree);  // ここではコスト引数は省略できる
```

`LinearGapCost(mismatch, gap)` — ギャップ 1 文字ごとに `gap` のコストがかかる。隣にギャップがあるかどうかは影響しない。置換と挿入・削除を同じ重みで扱いたくない場合に使う。gap を重くすると、クエリと長さの異なるエントリが不利になる。

```cpp
DTPatricia<DnaAlphabet, LinearGapCost> aligner(tree, LinearGapCost(1, 3));
```

`AffineGapCost(mismatch, gap_open, gap_extend)` — 連続する `L` 文字のギャップのコストは `gap_open + gap_extend * L` である。したがって、長い挿入・削除 1 箇所のほうが、短いものが複数あるより安くなる。挿入・削除がまとまって生じやすい生物配列の比較では、これが通常の選択肢である。

```cpp
DTPatricia<DnaAlphabet, AffineGapCost> aligner(tree, AffineGapCost(1, 2, 1));
```

`LinearGapCost` と `AffineGapCost` はパラメータを実行時に受け取るので、コストオブジェクトをアライナのコンストラクタに渡す必要がある。どちらも、正でないコストに対して `std::invalid_argument` を投げる。

どれを選んだ場合でも、返される `score` はそのモデルのもとでの厳密な最適コストである。

## 6. 3 種類のクエリ

3 つともクエリ文字列を受け取り、`std::vector<AlignmentResult>` を返す。

`ed_to_all(query)` — 辞書の**すべての**エントリまでの距離。本当に全件が必要なときに使う。枝刈りが一切効かないため、3 つの中で最も高コストである。

```cpp
auto all = aligner.ed_to_all("ACGT");
```

`ed_within_k(query, k)` — 距離が `k` 以下（境界値を含む）のエントリのみ。既定で選ぶべきはこれである。閾値が厳密な枝刈りを駆動するので、`k` が小さいほど探索は大幅に速くなる。

```cpp
auto near = aligner.ed_within_k("ACGT", 2);
```

`ed_pth_smallest(query, p)` — 最も近い `p` 件。探索はある距離のエントリをすべて見つけてからでないと打ち切ってよいか判断できないため、**打ち切り距離に同点があると結果が `p` 件を超えることがある**。これは意図的な仕様である。ちょうど `p` 件に切り詰めると、等しく良い候補の中から恣意的に選ぶことになってしまう。厳密に `p` 件が必要なら、呼び出し側で切り詰めること。

```cpp
auto top3 = aligner.ed_pth_smallest("ACGT", 3);
```

この 3 つでは表現できない停止規則——時間制限、閾値と件数の組み合わせ、それまでに集まった結果に対する条件など——が必要な場合は、`search_kernel` を独自の述語とともに直接呼ぶ。契約については [`search_kernel`](20_api_reference.ja.md#search_kernel) を、実例については [`examples/basic_example.cpp`](../examples/basic_example.cpp) の 6 節と 7 節を参照。

## 7. 知っておくべき入力の制約

以下は辞書とクエリの両方に当てはまる。

- 文字列に `'\0'` を含めてはならない。コード 0 は終端記号として予約されている。
- アルファベットに含まれない文字はすべて 1 つのキャッチオールコードに畳まれるため、それらの文字は互いに一致するものとして扱われる。`DnaAlphabet` のもとでは `N` と `X` は等しい。
- 大文字と小文字は区別されない。`DnaAlphabet` と `RyAlphabet` では `U` は `T` として扱われる。
- 辞書は空文字列や重複した文字列を含んでよい。重複はインデックスごとに別々の結果として報告される。

辞書を保持するのは木であり、アライナは木への**参照**を保持する。したがって木は、その上に構築したすべてのアライナより長く生存していなければならない。1 つの木を構築して、その上に（例えば異なるコストモデルで）複数のアライナを作るのは問題なく、コストも小さい。

## 8. 次に読むもの

- [`examples/basic_example.cpp`](../examples/basic_example.cpp) — すべてのアルファベットとコストの組み合わせを 1 ファイルにまとめた実行可能な例。`search_kernel` の独自述語も 2 つ含む。
- [API リファレンス](20_api_reference.ja.md) — 上記すべての正確な契約。
- [DT-Patricia の拡張](21_extending.ja.md) — 独自のアルファベット/コストポリシーの書き方。
