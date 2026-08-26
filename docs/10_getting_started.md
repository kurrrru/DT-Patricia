# Getting started

This page takes you from an empty project to a working DT-Patricia query. Read it top to
bottom; every snippet is complete enough to compile and run as-is.

For the exact signature and behaviour of anything mentioned here, see
[the API reference](20_api_reference.md). To plug in your own alphabet or cost model, see
[Extending DT-Patricia](21_extending.md).

## 1. Requirements

- A compiler with C++20 support. CI verifies GCC and Clang on Ubuntu.
- CMake 3.21 or later (3.25 or later if you want to use `CMakePresets.json`).

DT-Patricia is header-only and has no external library dependencies. Nothing is compiled into
a binary you have to link against — adding the headers to your include path is enough.

## 2. Adding it to your project

Pick whichever of these matches how you manage dependencies.

### Option A: vendor the sources (`add_subdirectory`)

Place the repository somewhere inside your project (for example `external/DT-Patricia`) and
add it from your `CMakeLists.txt`:

```cmake
add_subdirectory(external/DT-Patricia)
target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

Linking against the `dt_patricia::dt_patricia` target propagates both the include path and
the C++20 requirement, so you do not need to set `CMAKE_CXX_STANDARD` yourself.

When DT-Patricia is consumed this way it is not the top-level project, so its examples and
tests are not built.

### Option B: fetch it at configure time (`FetchContent`)

```cmake
include(FetchContent)
FetchContent_Declare(dt_patricia
    GIT_REPOSITORY https://github.com/kurrrru/DT-Patricia.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(dt_patricia)

target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

### Option C: install it, then `find_package`

Install the headers and the CMake package files once:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/prefix
cmake --install build
```

Then, from any other project:

```cmake
find_package(dt_patricia REQUIRED)
target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

If the install prefix is not one CMake searches by default, point it there when configuring
the consuming project:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/prefix
```

The installed package is versioned with `SameMajorVersion` compatibility, so
`find_package(dt_patricia 0.1 REQUIRED)` accepts any `0.x` install.

### Option D: no CMake at all

Add `include/` to your compiler's include path and compile with C++20:

```bash
g++ -std=c++20 -O2 -Iexternal/DT-Patricia/include my_app.cpp -o my_app
```

## 3. Your first query

There is a single entry-point header. Including it gives you everything in this page:

```cpp
#include <dt_patricia/dt_patricia.hpp>
```

Every use of the library follows the same three steps:

1. Build a `PatriciaTree` from your dictionary (a `std::vector<std::string>`).
2. Create a `DTPatricia` aligner over that tree.
3. Run a query against the aligner.

```cpp
#include <iostream>
#include <string>
#include <vector>

#include <dt_patricia/dt_patricia.hpp>

int main() {
    using namespace dt_patricia;

    std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG", "TGCA"};

    // 1) build the dictionary index
    PatriciaTree<DnaAlphabet> tree(targets);

    // 2) create an aligner over it
    DTPatricia<DnaAlphabet, UnitCost> aligner(tree);

    // 3) query: every dictionary entry within edit distance 1 of "ACGT"
    for (const AlignmentResult &r : aligner.ed_within_k("ACGT", 1)) {
        std::cout << r.score << "  " << targets[r.string_id] << "\n";
    }
}
```

Expected output (the three distance-1 hits may appear in any order among themselves):

```
0  ACGT
1  AAGT
1  ACG
1  ACGA
```

Two things about the result type are worth internalising now:

- `string_id` is an **index into the vector you passed to the tree**, in the order you passed
  it. The library never copies your strings back to you; you look them up yourself.
- `score` is the **exact** edit distance under the cost model you chose. It is not an
  estimate, and no candidate within the threshold is ever missed.

Results come back in non-decreasing order of `score`, and each dictionary entry appears at
most once.

### Building and running it

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/my_app
```

Build in `Release` (or at least with optimisation enabled). The search kernel is written to
be vectorised by the compiler, and an unoptimised build is dramatically slower.

## 4. Choosing an alphabet

The alphabet is a compile-time policy, and it decides **which characters are considered
distinct**. Two characters that map to the same code are treated as equal, so changing the
alphabet changes the meaning of the distance itself.

Three alphabets ship with the library:

| Alphabet | Distinguishes | Use it when |
| --- | --- | --- |
| `DnaAlphabet` | `A` `C` `G` `T`, plus one catch-all code | DNA/RNA sequences. `U` is folded into `T` |
| `RyAlphabet` | purine (`A`,`G`,`R`) vs pyrimidine (`C`,`T`,`U`,`Y`), plus a catch-all | you care only about the purine/pyrimidine pattern |
| `ProteinAlphabet` | the 20 standard amino acids, plus one catch-all | protein sequences |

Start with `DnaAlphabet` for nucleotides and `ProteinAlphabet` for peptides. Reach for
`RyAlphabet` only when collapsing to two classes is genuinely what you want — under it,
`GCAT` and `ACGT` are both `RYRY` and their distance is 0.

If none of these fits your data, writing one policy class is enough to support any alphabet;
see [Extending DT-Patricia](21_extending.md).

Note that the alphabet is part of the type: a `DTPatricia<Alphabet, Cost>` can only be built
over a `PatriciaTree<Alphabet>` with the *same* `Alphabet`. Mismatches are a compile error,
not a silent wrong answer.

## 5. Choosing a cost model

The cost model decides **how much each edit operation costs**. Three are provided.

`UnitCost` — mismatch and gap both cost 1, i.e. plain Levenshtein distance. This is the
default, and it is also the fastest: it takes a specialised code path that the other two
cannot use. Use it unless you have a specific reason not to.

```cpp
DTPatricia<DnaAlphabet, UnitCost> aligner(tree);  // the cost argument is optional here
```

`LinearGapCost(mismatch, gap)` — each gap character costs `gap`, independently of whether it
is adjacent to another gap. Use it when substitutions and indels should not weigh the same.
Making the gap heavier penalises entries whose length differs from the query.

```cpp
DTPatricia<DnaAlphabet, LinearGapCost> aligner(tree, LinearGapCost(1, 3));
```

`AffineGapCost(mismatch, gap_open, gap_extend)` — a run of `L` consecutive gap characters
costs `gap_open + gap_extend * L`. One long indel is therefore cheaper than several short
ones. This is the usual choice in biological sequence comparison, where indels tend to occur
as runs.

```cpp
DTPatricia<DnaAlphabet, AffineGapCost> aligner(tree, AffineGapCost(1, 2, 1));
```

`LinearGapCost` and `AffineGapCost` take their parameters at run time, so the cost object has
to be passed to the aligner's constructor. Both reject non-positive costs by throwing
`std::invalid_argument`.

Whichever you pick, the returned `score` is the exact optimal cost under that model.

## 6. The three queries

All three take the query string and return `std::vector<AlignmentResult>`.

`ed_to_all(query)` — the distance to **every** entry in the dictionary. Use it when you
really do want all of them; it is the most expensive of the three because nothing can be
pruned.

```cpp
auto all = aligner.ed_to_all("ACGT");
```

`ed_within_k(query, k)` — only the entries whose distance is at most `k`, boundary included.
This is the one to reach for by default: the threshold drives an exact pruning step, so a
small `k` makes the search much faster.

```cpp
auto near = aligner.ed_within_k("ACGT", 2);
```

`ed_kth_smallest(query, k)` — the `k` nearest entries. Because the search finds every entry
at a given distance before it can know whether it has enough, **the result may contain more
than `k` entries when there are ties at the cut-off distance**. That is deliberate:
truncating to exactly `k` would mean picking arbitrarily among entries that are equally
good. Truncate yourself if your application needs exactly `k`.

```cpp
auto top3 = aligner.ed_kth_smallest("ACGT", 3);
```

If you need a stopping rule these three do not express — a time budget, a combination of a
threshold and a count, a condition on the results collected so far — call `search_kernel`
directly with your own predicate. See
[`search_kernel`](20_api_reference.md#search_kernel) for the contract, and sections 6 and 7
of [`examples/basic_example.cpp`](../examples/basic_example.cpp) for worked examples.

## 7. Input constraints to know about

These apply to both the dictionary and the query:

- Strings must not contain `'\0'`. Code 0 is reserved as the terminator.
- Every character outside the alphabet collapses onto a single catch-all code, so all such
  characters are treated as matching one another. Under `DnaAlphabet`, `N` and `X` are equal.
- Case is not distinguished. In `DnaAlphabet` and `RyAlphabet`, `U` is treated as `T`.
- The dictionary may contain empty strings and duplicate strings. Duplicates are reported as
  separate results, one per index.

The tree holds the dictionary; the aligner holds a **reference** to the tree. The tree must
therefore outlive every aligner built over it. Building one tree and creating several
aligners over it (with different cost models, say) is fine and cheap.

## 8. Where to go next

- [`examples/basic_example.cpp`](../examples/basic_example.cpp) — every alphabet and cost
  combination in one runnable file, plus two custom `search_kernel` predicates.
- [API reference](20_api_reference.md) — the exact contract of everything above.
- [Extending DT-Patricia](21_extending.md) — writing your own alphabet or cost policy.
