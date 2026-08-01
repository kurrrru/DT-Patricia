[![tests](https://github.com/kurrrru/DT-Patricia/actions/workflows/tests.yml/badge.svg)](https://github.com/kurrrru/DT-Patricia/actions/workflows/tests.yml)
[![format](https://github.com/kurrrru/DT-Patricia/actions/workflows/format.yml/badge.svg)](https://github.com/kurrrru/DT-Patricia/actions/workflows/format.yml)

> [!WARNING]
> **This repository is a work in progress.** The interfaces of its functions and classes, the internal implementation, and the directory layout may change without notice.

# DT-Patricia

**DT-Patricia** (*Diagonal Transition on Patricia Tree*) is a header-only C++20 library that computes the exact edit distance between a single query sequence and every string in a dictionary (a set of strings).

The dictionary is stored in a Patricia tree (a path-compressed trie), and **the diagonal transition algorithm is run on that tree**. Instead of processing the dictionary one string at a time, the basic idea is to perform the computation for a shared prefix once.

## Features

- **Exact results** — the returned score is always the optimal edit distance; no approximations or heuristics are involved
- **Three cost models** — unit cost / linear gap cost / affine gap cost
- **Pluggable alphabets** — DNA (A/C/G/T), R/Y (purine / pyrimidine), and protein (20 amino acids) are bundled. Writing a single policy class is enough to support any alphabet
- **Three kinds of queries** — the score for every entry / entries within distance k / the k nearest entries
- **Header-only, no external dependencies** — just include it on any toolchain with C++20

## Requirements

- A compiler with C++20 support (CI verifies GCC and Clang on Ubuntu)
- CMake 3.21 or later (3.25 or later when using `CMakePresets.json`)

There are no external library dependencies.

## Quick start

```cpp
#include <iostream>
#include <string>
#include <vector>

#include <dt_patricia/dt_patricia.hpp>

int main() {
    using namespace dt_patricia;

    std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG", "TGCA"};

    PatriciaTree<DnaAlphabet> tree(targets);          // 1) build a Patricia tree from the dictionary
    DTPatricia<DnaAlphabet, UnitCost> aligner(tree);  // 2) create an aligner

    for (const auto &r : aligner.ed_within_k("ACGT", 1)) {  // 3) query
        std::cout << r.score << "  " << targets[r.string_id] << "\n";
    }
}
```

For examples that combine alphabets with cost models, and for searches driven by a custom stop condition, see [`examples/basic_example.cpp`](examples/basic_example.cpp).

### Building and running the example

```bash
# 1. configure: -S = where the sources are (.), -B = where to build (build)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# 2. build
cmake --build build -j

# 3. run
./build/examples/basic_example
```

Or use the presets. Two are available: `default` (Release, `build/default/`) and `debug` (Debug, `build/debug/`).

```bash
cmake --preset default
cmake --build --preset default -j
./build/default/examples/basic_example
```

### Using it in your project

Linking against the CMake target `dt_patricia::dt_patricia` propagates both the include path and the C++20 requirement.

```cmake
add_subdirectory(external/DT-Patricia)
target_link_libraries(my_app PRIVATE dt_patricia::dt_patricia)
```

If you do not use CMake, adding `include/` to your include path is enough.

### Tests

The tests check that the results agree exactly with a brute-force implementation based on naive DP.

Start by generating the random test cases (`tests/testcase_random/` is not tracked by Git, so it has to be generated before running them).

```bash
# defaults: seed 42, 60 cases, written to tests/testcase_random
python3 tests/scripts/gen_random_tests.py

# the seed, the number of cases, and the output directory can be changed
python3 tests/scripts/gen_random_tests.py --seed 7 --count 300
python3 tests/scripts/gen_random_tests.py --out-dir /path/to/dir

# by default the output directory is emptied first; pass --no-clean to add to it
python3 tests/scripts/gen_random_tests.py --no-clean
```

The generated files are fully reproducible for the same seed and count.

Build and run them as follows.

```bash
cmake -S . -B build/tests \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DDT_PATRICIA_BUILD_TESTS=ON \
  -DDT_PATRICIA_BUILD_EXAMPLES=OFF
cmake --build build/tests -j
ctest --test-dir build/tests --output-on-failure
```

Or use the preset.

```bash
cmake --preset tests             # configure -> generated into build/tests/
cmake --build --preset tests -j  # build
ctest --preset tests -j          # run (--output-on-failure is built into the preset)
```

## Input constraints

- Strings must not contain `'\0'` (code 0 is reserved as the terminator).
- Every character outside the alphabet collapses onto a single code, so all such characters are treated as matching one another.
- Case is not distinguished. In `DnaAlphabet` and `RyAlphabet`, the RNA character `'U'` is treated as `'T'`.
- The dictionary may contain empty strings and duplicate strings. Duplicates are returned as separate results.

## Directory layout

```
DT-Patricia/
├── include/dt_patricia/            # the library itself (all headers)
│   ├── dt_patricia.hpp             # the only entry point
│   ├── patricia_tree.hpp           # the Patricia tree that holds the dictionary
│   ├── aligner.hpp                 # the search engine, DTPatricia
│   ├── alignment_result.hpp        # the result type, AlignmentResult
│   ├── policy/                     # alphabets and cost models
│   ├── internal/                   # internal implementation
│   └── debug/                      # debugging helpers
├── examples/                       # usage examples
├── tests/                          # test driver, test cases, generator script
├── cmake/                          # CMake package files used on install
└── external/                       # place for external dependencies (currently empty)
```
