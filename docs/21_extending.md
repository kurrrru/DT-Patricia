# Extending DT-Patricia

DT-Patricia is parameterised over two policies: an **alphabet**, which decides which
characters are considered distinct, and a **cost model**, which decides what each edit
operation costs. Both are plain structs you can write yourself, and neither requires touching
the library's internals.

This page is the normative specification of what those structs must provide. For advice on
which of the bundled policies to pick, see [Getting started](10_getting_started.md); for the
behaviour of the bundled ones, see the [API reference](20_api_reference.md).

## `AlphabetPolicy`

```cpp
#include <dt_patricia/policy/alphabet.hpp>

template <class Alphabet>
concept AlphabetPolicy = /* see below */;
```

An alphabet policy defines a mapping from all 256 byte values onto a compact range of codes
`0 .. CODE_MAX`. The library never compares characters directly: it compares codes. Two
characters mapped to the same code are therefore indistinguishable — that is the mechanism
behind `RyAlphabet` folding `A` and `G` together, and behind every alphabet's catch-all code.

Both the tree and the aligner are templated on the policy, so all of this is resolved at
compile time and the tables are baked into the binary.

### Required members

A policy is a class or struct with exactly four public static members. It holds no state and
is never instantiated.

```cpp
static constexpr uint8_t CODE_TERM;
static constexpr uint8_t CODE_MAX;
static consteval std::array<uint8_t, 256> make_char_to_code() noexcept;
static consteval std::array<char, CODE_MAX + 1> make_code_to_char() noexcept;
```

- **`CODE_TERM`** — The terminator code. **Must be 0.** It is fixed rather than configurable
  because the tree's node layout relies on the terminator occupying the first bucket.
- **`CODE_MAX`** — The largest code the policy ever produces. **Must be greater than 0.** The
  tree allocates `CODE_MAX + 1` transition slots per node, so this value directly determines
  memory use: keep the code range tight and contiguous rather than, say, using ASCII values as
  codes.
- **`make_char_to_code()`** — Returns the encoding table, indexed by `unsigned char`. Every one
  of the 256 entries must be filled, and every value must lie in `0 .. CODE_MAX`. This is the
  table that folds case, synonyms and unknown characters together.
- **`make_code_to_char()`** — Returns the inverse table: for each code, the single character
  that represents it. Used to canonicalise strings and to render labels.

Both table functions are `consteval` and `noexcept`, and both must return exactly the types
above — the concept checks the return type with `std::same_as`, so `auto` deduction that
yields a different array size will not satisfy it.

### Invariants

The concept does not stop at checking that the members exist. It also evaluates
`detail::validate_alphabet_policy<Alphabet>()` at compile time, which enforces:

1. `CODE_TERM == 0`.
2. `CODE_MAX > 0`.
3. `make_char_to_code()['\0'] == CODE_TERM`.
4. Every entry of `make_char_to_code()` is `<= CODE_MAX`.
5. `make_code_to_char()[CODE_TERM] == '\0'`.
6. **Round-trip consistency:** for every code `c` in `0 .. CODE_MAX`,
   `make_char_to_code()[make_code_to_char()[c]] == c`.

Rule 6 is the one that catches real mistakes. It says every code must have a representative
character that encodes back to that same code — which in particular means no code may be
unused, and no two codes may share a representative. If you add a code to one table and
forget the other, the policy stops satisfying the concept and instantiating
`PatriciaTree` with it fails to compile.

Because the whole check runs at compile time, a broken policy can never produce wrong
answers at run time; it simply does not build.

### Design rules that the concept cannot check

- **Give unknown characters a catch-all code.** The tables must cover all 256 byte values, so
  every character you do not name has to go *somewhere*. Reserve one code for them, as the
  bundled policies do. Be aware of the consequence: all such characters then match one
  another. Under `ProteinAlphabet`, `B` and `Z` are equal, because both are catch-all.
- **Never map anything but `'\0'` to code 0.** Code 0 is the terminator; a real character
  sharing it would corrupt the tree.
- **Do not put per-character code constants in the public interface.** The bundled policies
  keep `CODE_A`, `CODE_C` and friends private, and the concept deliberately does not require
  them. Keeping them private means the numeric values remain free to change.
- **Keep codes dense.** A gap in `0 .. CODE_MAX` is rejected by rule 6 anyway, but the wider
  point is that every extra code costs a transition slot on every node.

### Worked example

An alphabet for ordinary lowercase text — 26 letters, case-insensitive, everything else
collapsed:

```cpp
#include <array>
#include <cstdint>

#include <dt_patricia/dt_patricia.hpp>

struct LowercaseAlphabet {
 public:
    static constexpr uint8_t CODE_TERM = 0;
    static constexpr uint8_t CODE_MAX = 27;  // 0 = terminator, 1..26 = 'a'..'z', 27 = other

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
        table[CODE_OTHER] = '?';  // the representative of everything else
        return table;
    }

 private:
    static constexpr uint8_t CODE_OTHER = 27;
};

static_assert(dt_patricia::AlphabetPolicy<LowercaseAlphabet>);
```

The `static_assert` is not required, but put one next to every policy you write: it reports
the violation at the policy's own definition rather than deep inside a `PatriciaTree`
instantiation.

Using it is no different from using a bundled policy:

```cpp
std::vector<std::string> dictionary = {"receive", "recieve", "retrieve"};
dt_patricia::PatriciaTree<LowercaseAlphabet> tree(dictionary);
dt_patricia::DTPatricia<LowercaseAlphabet, dt_patricia::UnitCost> aligner(tree);

auto hits = aligner.ed_within_k("recieve", 2);  // finds "receive" too
```

## Cost policies

A cost policy is deliberately **not** constrained by a named concept: it is duck-typed, and
which members are required depends on the flags the policy itself declares.

Unlike an alphabet policy, a cost policy carries run-time state — its cost values are
ordinary non-static members — and the aligner stores a copy of it.

### Required members

Every cost policy must provide two compile-time flags:

```cpp
static constexpr bool is_linear;
static constexpr bool is_unit;
```

- **`is_linear`** — `true` if a run of `L` gap characters costs `gap * L`. `false` selects the
  affine model, where opening a run costs extra. This flag chooses between the two
  `search_kernel` overloads, which are separate implementations.
- **`is_unit`** — `true` selects a specialised fast path. **Only set it when the model is
  exactly Levenshtein**, i.e. `is_linear == true` and both `mismatch` and `gap` equal 1. The
  fast path does not read the cost values at all, so declaring `is_unit = true` alongside any
  other values silently computes unit-cost distances instead of yours. This is the one
  requirement here that nothing checks for you.

Beyond the flags, the members required depend on `is_linear`:

| `is_linear` | Required members | Cost of a gap run of length `L` |
| --- | --- | --- |
| `true` | `mismatch`, `gap` | `gap * L` |
| `false` | `mismatch`, `gap_open`, `gap_extend` | `gap_open + gap_extend * L` |

All cost members are read as `uint32_t`. They may be `static constexpr` (as in `UnitCost`) or
ordinary members set at construction (as in `LinearGapCost` and `AffineGapCost`); the
implementation only ever reads them through the stored cost object, so either works.

### Invariants

These are your responsibility to uphold. The bundled policies enforce them in their
constructors by throwing `std::invalid_argument`; a policy of your own should do the same.

- `mismatch >= 1`.
- `gap >= 1` when `is_linear`. The pruning step divides by `gap`, so 0 is not merely
  meaningless but fatal.
- `gap_extend >= 1` when `!is_linear`. `gap_open` may be 0, which degenerates to the linear
  model.
- `is_unit` implies `mismatch == gap == 1`.

Costs are also what determines how much wavefront history the search keeps: the number of
retained levels grows with the largest cost value. Very large costs are correct but consume
proportionally more memory, so prefer the smallest integers that express the ratio you want —
`LinearGapCost(1, 3)` rather than `LinearGapCost(100, 300)`.

### Worked example

A linear model with a fixed, compile-time cost ratio. Because everything is `static
constexpr`, the type is default-constructible and can be used without passing a cost object:

```cpp
struct GapHeavyCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit = false;  // NOT unit: the costs are not 1/1
    static constexpr uint32_t mismatch = 1;
    static constexpr uint32_t gap = 3;
};

dt_patricia::DTPatricia<dt_patricia::DnaAlphabet, GapHeavyCost> aligner(tree);
```

A run-time-configurable one follows the shape of `LinearGapCost`:

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

Since this type has no default constructor, the cost argument is mandatory — the aligner's
default argument is only instantiated when you omit it.

## Verifying a new policy

The library's own tests compare its results against a brute-force dynamic-programming
implementation, and that is the right way to gain confidence in a policy of your own:
generate random dictionaries and queries over your alphabet, compute distances with a naive
`O(nm)` DP using the same cost values, and check the two agree exactly. See `tests/` for the
existing harness.

For an alphabet policy specifically, `canonicalize<Alphabet>` is the quickest sanity check:
it shows exactly what the library will compare, so `canonicalize<LowercaseAlphabet>("Hello!")`
returning `"hello?"` confirms the folding rules in one line.
