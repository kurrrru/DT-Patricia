# API reference

The exact contract of DT-Patricia's public interface. If you are looking for how to
accomplish something rather than what a given entity does, start at
[Getting started](10_getting_started.md).

Everything below lives in namespace `dt_patricia`.

## Contents

- [Headers](#headers)
- [`AlignmentResult`](#alignmentresult)
- [`PatriciaTree<Alphabet>`](#patriciatreealphabet)
- [`DTPatricia<Alphabet, CostType>`](#dtpatriciaalphabet-costtype)
- [Query result semantics](#query-result-semantics)
- [Cost policies](#cost-policies)
- [Alphabet policies](#alphabet-policies)
- [Free functions](#free-functions)
- [Thread safety](#thread-safety)
- [Exceptions](#exceptions)

## Headers

| Header | Contents |
| --- | --- |
| `<dt_patricia/dt_patricia.hpp>` | Entry point. Includes everything below. Include this one. |
| `<dt_patricia/patricia_tree.hpp>` | `PatriciaTree` |
| `<dt_patricia/aligner.hpp>` | `DTPatricia` |
| `<dt_patricia/alignment_result.hpp>` | `AlignmentResult` |
| `<dt_patricia/policy/alphabet.hpp>` | `AlphabetPolicy`, the bundled alphabets, `canonicalize` |
| `<dt_patricia/policy/cost.hpp>` | `UnitCost`, `LinearGapCost`, `AffineGapCost` |

`dt_patricia/internal/` and `dt_patricia/debug/` are not part of the public interface. Their
contents may change at any time, so including them or naming anything in
`dt_patricia::internal` is not recommended.

## `AlignmentResult`

```cpp
struct AlignmentResult {
    uint32_t string_id;
    uint32_t score;
};
```

- **`string_id`** — Index of the matched string in the `std::vector<std::string>` that was
  passed to the `PatriciaTree` constructor, in the original input order. Not an internal node
  id.
- **`score`** — The exact optimal alignment cost between the query and that string, under the
  aligner's cost model. With `UnitCost` this is the Levenshtein distance.

An aggregate; it has no member functions.

## `PatriciaTree<Alphabet>`

```cpp
template <AlphabetPolicy Alphabet = DnaAlphabet>
class PatriciaTree;
```

An immutable, path-compressed trie holding the dictionary. Build one, then create one or
more aligners over it.

`Alphabet` must satisfy the [`AlphabetPolicy`](21_extending.md#alphabetpolicy) concept.

### Construction and special members

```cpp
explicit PatriciaTree(const std::vector<std::string> &input_data);

PatriciaTree() = delete;
PatriciaTree(const PatriciaTree &) = delete;             // dictionaries can be large
PatriciaTree &operator=(const PatriciaTree &) = delete;
PatriciaTree(PatriciaTree &&) noexcept = default;
PatriciaTree &operator=(PatriciaTree &&) noexcept = default;
~PatriciaTree() = default;
```

The constructor builds the whole tree; there is no separate build step and no way to insert
or erase afterwards. `input_data` is not retained — the tree copies what it needs — so it may
be destroyed as soon as the constructor returns. String ids refer to positions in
`input_data`, so callers who want the strings back must keep them.

`input_data` may be empty, may contain empty strings, and may contain duplicates. Duplicates
keep their distinct ids and are reported as separate results.

No string may contain `'\0'`.

The tree is copyable-by-construction only: it is movable but not copyable.

### Member types and constants

```cpp
using alphabet_type = Alphabet;

static constexpr uint8_t     CODE_TERM = Alphabet::CODE_TERM;   // always 0
static constexpr uint8_t     CODE_MAX  = Alphabet::CODE_MAX;
static constexpr std::size_t BUCKET_SIZE = CODE_MAX + 1;
static constexpr std::size_t SIMD_PADDING_SIZE = 256;

inline static constexpr std::array<uint8_t, 256> CHAR_TO_CODE = Alphabet::make_char_to_code();
```

`SIMD_PADDING_SIZE` is the slack the search kernel appends to its internal copy of the query
so that vectorised loads cannot read past the end. It is exposed as a constant rather than as
a requirement on the caller: query strings you pass in need no padding of their own.

### Dictionary-level queries

```cpp
[[nodiscard]] bool     empty() const noexcept;
[[nodiscard]] uint32_t string_count() const noexcept;
[[nodiscard]] uint32_t node_count() const noexcept;
```

- **`empty()`** — `true` if the tree holds no nodes, i.e. it was built from an empty vector.
  Querying an empty tree is well-defined and returns no results.
- **`string_count()`** — Number of strings in the dictionary, counting duplicates separately.
  Equal to `input_data.size()`.
- **`node_count()`** — Number of node slots in the tree, including the reserved id 0. An
  implementation detail of the layout rather than a meaningful measure of dictionary size.

### Node-level traversal

These expose the tree structure itself. They are public because they are useful for
inspecting, exporting or reimplementing traversals over the dictionary, but a normal user of
the library never needs them — the search in `DTPatricia` is the intended way to consume the
tree.

Node id `0` is the reserved invalid value; the root is id `1`. Every accessor below expects a
valid node id; passing `0` or an out-of-range id is undefined behaviour except where noted.

```cpp
[[nodiscard]] uint32_t root_id() const noexcept;                              // always 1
[[nodiscard]] uint32_t transition(uint32_t node_id, char ch) const noexcept;
[[nodiscard]] uint32_t transition(uint32_t node_id, uint8_t code) const noexcept;
[[nodiscard]] uint32_t get_parent(uint32_t node_id) const noexcept;
[[nodiscard]] bool     is_leaf(uint32_t node_id) const noexcept;
[[nodiscard]] bool     is_terminal(uint32_t node_id) const noexcept;
[[nodiscard]] std::string_view get_label(uint32_t node_id) const noexcept;
[[nodiscard]] uint32_t get_label_length(uint32_t node_id) const noexcept;
[[nodiscard]] std::span<const uint32_t> get_string_id(uint32_t node_id) const noexcept;
```

- **`transition(node_id, ch)`** — Follows the edge leaving `node_id` whose first character is
  `ch`, and returns the id of the node it reaches, or `0` if there is no such edge. The `char`
  overload encodes `ch` through `CHAR_TO_CODE` first, so it applies the alphabet's folding rules
  (case, `U`/`T`, catch-all). The `uint8_t` overload takes an already-encoded code.
- **`get_parent(node_id)`** — The parent's node id, or `0` for the root and for out-of-range
  ids.
- **`is_leaf(node_id)`** — `true` if the node has no outgoing edges.
- **`is_terminal(node_id)`** — `true` if a dictionary string ends at this node, i.e. it has an
  outgoing `CODE_TERM` edge. Note that the string ids are **not** attached to the terminal node
  itself but to the leaf reached through that edge — `get_string_id(transition(node_id,
  CODE_TERM))`.
- **`get_label(node_id)`** — The label of the edge **entering** `node_id`. The returned view
  points into storage owned by the tree and is valid as long as the tree is alive and not moved
  from.
- **`get_string_id(node_id)`** — The ids of the dictionary strings ending at this leaf. Returns
  an empty span if `node_id` is not a leaf. More than one id is returned when the dictionary
  contained duplicates.

### Precomputed subtree statistics

```cpp
[[nodiscard]] const std::vector<uint32_t> &get_subtree_counts() const noexcept;
[[nodiscard]] const std::vector<uint32_t> &get_subtree_max_lengths() const noexcept;
[[nodiscard]] const std::vector<uint32_t> &get_subtree_min_lengths() const noexcept;
[[nodiscard]] const std::vector<uint32_t> &get_parent_path_lengths() const noexcept;
```

Each is indexed by node id and is computed once during construction. The search uses them to
prune; they are exposed for the same reason as the traversal accessors above.

- **`get_subtree_counts()`** — Number of dictionary strings in the subtree rooted at each node.
- **`get_subtree_max_lengths()` / `get_subtree_min_lengths()`** — The longest and shortest
  distance from each node down to a leaf of its subtree, counting that node's own incoming
  label.
- **`get_parent_path_lengths()`** — The length of the path from the root to each node's
  **parent**, excluding the node's own incoming label.

The returned references are valid as long as the tree is alive and not moved from.

## `DTPatricia<Alphabet, CostType>`

```cpp
template <AlphabetPolicy Alphabet = DnaAlphabet, typename CostType = UnitCost>
class DTPatricia;
```

The search engine. It runs the diagonal transition algorithm over a `PatriciaTree`.

`Alphabet` must be the *same* type the tree was instantiated with; a mismatch is a compile
error. `CostType` must satisfy the [cost policy requirements](21_extending.md#cost-policies).

### Construction and special members

```cpp
DTPatricia(const tree_type &patricia_tree, CostType cost = CostType());

DTPatricia() = delete;
DTPatricia(const DTPatricia &) = delete;
DTPatricia &operator=(const DTPatricia &) = delete;
DTPatricia(DTPatricia &&) noexcept = delete;
DTPatricia &operator=(DTPatricia &&) noexcept = delete;
~DTPatricia() = default;
```

The aligner stores a **reference** to the tree and a **copy** of the cost object. The tree
must outlive the aligner. Construction is cheap: no per-aligner index is built, so it is
possible to create several aligners with different cost models over one tree.

Neither copyable nor movable, precisely because it holds that reference.

The default argument for `cost` is only instantiated if you omit the argument, so a
`CostType` without a default constructor — `LinearGapCost` and `AffineGapCost` are both such
types — simply requires you to pass one.

### Member types

```cpp
using alphabet_type = Alphabet;
using tree_type     = PatriciaTree<Alphabet>;
```

### `get_patricia_tree`

```cpp
[[nodiscard]] const tree_type &get_patricia_tree() const noexcept;
```

The tree the aligner was constructed over.

### `ed_to_all`

```cpp
std::vector<AlignmentResult> ed_to_all(const std::string &query) const;
```

Returns the distance from `query` to **every** string in the dictionary — exactly
`string_count()` results, or none if the tree is empty.

No threshold is available to prune with, so this is the most expensive of the three queries.

### `ed_within_k`

```cpp
std::vector<AlignmentResult> ed_within_k(const std::string &query, int k) const;
```

Returns every dictionary string whose distance from `query` is **at most `k`**, boundary
included. `k` is a cost under the aligner's cost model, not a number of edit operations, so
under `LinearGapCost(1, 3)` a single gap already costs 3.

`k` doubles as an upper bound for pruning, and the pruning is exact with respect to it: no
string within the threshold is ever dropped. Smaller `k` therefore means a faster search.

Because a distance level is always completed before the stop condition is evaluated, `k == 0` returns the exact matches, and any `k < 0` behaves the same as `k == 0` rather than returning nothing.

### `ed_kth_smallest`

```cpp
std::vector<AlignmentResult> ed_kth_smallest(const std::string &query, size_t k) const;
```

Returns the `k` nearest dictionary strings. `k` is clamped to `string_count()`. For the same reason as above, `k == 0` returns the exact matches rather than nothing.

**The result may contain more than `k` entries.** The search completes each distance level
before it can test whether enough results have been collected, so if several strings tie at
the cut-off distance, all of them are returned. This avoids choosing arbitrarily between
equally good candidates. Truncate the result yourself if you need exactly `k`.

No pruning bound is used, because the cut-off distance is not known in advance.

### `search_kernel`

```cpp
template <typename StopPredicate>
std::vector<AlignmentResult> search_kernel(const std::string &query,
                                           StopPredicate stop_predicate,
                                           int upper_bound = -1) const;
```

The engine the other three queries are thin wrappers around. Use it when you need a stopping
rule they do not express.

- **`stop_predicate`** — Callable with the signature `bool(int current_score, const
  std::vector<AlignmentResult> &results)`. It is invoked **once per distance level**, after
  every match at that distance has been appended to `results`. Returning `true` stops the search
  and the collected results are returned as-is. `current_score` is the distance level just
  finished, and increases monotonically across calls. Because it is consulted only between
  levels, any budget expressed through it — including a wall-clock deadline — has a granularity
  of one distance level's worth of work.
- **`upper_bound`** — Enables pruning at cost `upper_bound`, exactly as `ed_within_k` does.
  `-1`, the default, disables pruning. The pruning never discards a candidate that could be
  within the bound, but note that it does **not** by itself stop the search — pass a
  `stop_predicate` that also respects the bound, or the search will keep going past it with
  nothing left to find.

Two overloads exist, constrained on `CostType::is_linear`; the correct one is selected
automatically and they are indistinguishable to the caller.

Sections 6 and 7 of [`examples/basic_example.cpp`](../examples/basic_example.cpp) show a
timeout predicate and a combined top-k-with-threshold predicate.

## Query result semantics

These hold for all four query functions.

- **Order.** Results are appended in non-decreasing order of `score`. The relative order of
  entries sharing a score is unspecified and should not be relied upon.
- **Uniqueness.** Each `string_id` appears at most once. Duplicate strings in the dictionary
  have distinct ids, so each of them appears once.
- **Exactness.** Every reported `score` is the optimal cost. No approximation, no heuristic.
- **Empty tree.** Returns an empty vector.
- **Empty query.** Well-defined: the distance to each dictionary string is the cost of
  deleting it entirely.
- **Query length.** The query is copied and canonicalised internally; the string you pass is
  not modified.

## Cost policies

Defined in `<dt_patricia/policy/cost.hpp>`. All costs are `uint32_t`.

### `UnitCost`

```cpp
struct UnitCost {
    static constexpr bool is_linear = true;
    static constexpr bool is_unit   = true;
    static constexpr uint32_t mismatch = 1;
    static constexpr uint32_t gap      = 1;
};
```

Plain Levenshtein distance. Stateless and default-constructible, so
`DTPatricia<Alphabet, UnitCost>` needs no cost argument. `is_unit` selects a specialised code
path, making this materially faster than expressing the same costs through `LinearGapCost`.

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

A run of `L` gap characters costs `gap * L`. Throws `std::invalid_argument` if `m < 1` or
`g < 1`. Not default-constructible.

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

A run of `L` gap characters costs `gap_open + gap_extend * L`. Throws
`std::invalid_argument` if `m < 1` or `g_extend < 1`; `gap_open` may be 0. Not
default-constructible.

## Alphabet policies

Defined in `<dt_patricia/policy/alphabet.hpp>`. Each maps all 256 byte values onto a small
set of codes. Code `0` is always the terminator. Characters sharing a code are
indistinguishable to the search, and each alphabet has one catch-all code that every
character it does not name falls into — so all unnamed characters match one another.

Case is never distinguished.

| Policy | `CODE_MAX` | Named characters | Catch-all representative |
| --- | --- | --- | --- |
| `DnaAlphabet` | 5 | `A` `C` `G` `T`; `U` folds to `T` | `N` |
| `RyAlphabet` | 3 | `R` (from `A`, `G`, `R`), `Y` (from `C`, `T`, `U`, `Y`) | `N` |
| `ProteinAlphabet` | 21 | `A C D E F G H I K L M N P Q R S T V W Y` | `X` |

Under `ProteinAlphabet` the ambiguity codes `B`, `J`, `O`, `U`, `X` and `Z` all fall into the
catch-all and are therefore mutually equal.

The `AlphabetPolicy` concept these satisfy is documented in
[Extending DT-Patricia](21_extending.md#alphabetpolicy).

## Free functions

```cpp
template <AlphabetPolicy Alphabet>
std::string canonicalize(std::string_view s);

template <AlphabetPolicy Alphabet>
void canonicalize_inplace(char *data, std::size_t n) noexcept;
```

Rewrite each character to the representative character of its code under `Alphabet`. This is
what the search does internally to the query, and it is exposed so that callers can see
exactly what the library compares — under `DnaAlphabet`, `canonicalize<DnaAlphabet>("acgu")`
is `"ACGT"`, and under `RyAlphabet` it is `"RYRY"`.

Two strings have distance 0 if and only if their canonical forms are equal.

`canonicalize_inplace` overwrites `n` bytes starting at `data`. The lookup table is a
compile-time constant with static storage duration, so both functions are safe to call from
multiple threads.

## Thread safety

`PatriciaTree` is immutable after construction, and every query function on `DTPatricia` is
`const` and keeps all of its working state in local variables. Any number of threads may
therefore query the same tree and even the same aligner concurrently, with no external
synchronisation.

There is no internal parallelism: one query runs on one thread.

## Exceptions

The library throws only from the cost policy constructors, which raise
`std::invalid_argument` on non-positive costs. Allocation failures from the standard library
propagate normally. The query functions themselves do not throw.
