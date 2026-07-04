// A sample that demonstrates the basic usage of DT-Patricia.
// By switching the "alphabet (DNA / R-Y / Protein)" and the "cost
// (Unit / Linear / Affine)", you can see how the results change for the
// same query and dictionary.
//
//   The usage is always the same 3 steps:
//     1) Build a PatriciaTree<Alphabet> from a std::vector<std::string>
//     2) Pass the tree (and a cost) to DTPatricia<Alphabet, Cost>
//     3) Query with ed_to_all / ed_within_k / ed_kth_smallest,
//        or call search_kernel directly with your own stop predicate

#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <dt_patricia/dt_patricia.hpp>

namespace {

void print_header(const std::string &title) {
    std::cout << "\n========== " << title << " ==========" << std::endl;
}

// Print the matches. The fields of AlignmentResult are:
//   string_id : index into the original `targets` dictionary (input order)
//   score     : edit distance (cost) from the query
// The results come back sorted by ascending score, so we print them as-is.
void print_results(const std::vector<std::string> &targets, const std::string &query,
                   const std::vector<dt_patricia::AlignmentResult> &results) {
    std::cout << "  query = \"" << query << "\"  (" << results.size() << " hit)" << std::endl;
    for (const auto &r : results) {
        std::cout << "    score=" << std::setw(2) << r.score << "  id=" << std::setw(2)
                  << r.string_id << "  target=\"" << targets[r.string_id] << "\"" << std::endl;
    }
}

}  // namespace

int main() {
    using namespace dt_patricia;

    // ============================================================
    // 1) DNA + UnitCost
    //    The most basic case. Both mismatch and gap cost 1, i.e. the
    //    plain edit distance. ed_to_all returns the distance to every
    //    entry in the dictionary.
    // ============================================================
    {
        print_header("DNA / UnitCost / ed_to_all");
        std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG", "TGCA", "GGGG"};
        PatriciaTree<DnaAlphabet> tree(targets);
        DTPatricia<DnaAlphabet, UnitCost> aligner(tree);  // 2nd arg (cost) is optional

        const std::string query = "ACGT";
        print_results(targets, query, aligner.ed_to_all(query));
    }

    // ============================================================
    // 2) DNA + LinearGapCost(mismatch=1, gap=3)
    //    Making the gap heavier raises the score of entries that require
    //    a length change (insertion/deletion). Compared with UnitCost,
    //    entries such as "ACG" change.
    // ============================================================
    {
        print_header("DNA / LinearGapCost(mismatch=1, gap=3) / ed_to_all");
        std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG", "TGCA", "GGGG"};
        PatriciaTree<DnaAlphabet> tree(targets);
        DTPatricia<DnaAlphabet, LinearGapCost> aligner(tree, LinearGapCost(1, 3));

        const std::string query = "ACGT";
        print_results(targets, query, aligner.ed_to_all(query));
    }

    // ============================================================
    // 3) DNA + AffineGapCost(mismatch=1, gap_open=2, gap_extend=1)
    //    Gap cost = gap_open + gap_extend * L. Opening consecutive gaps as
    //    a single run is cheaper than opening them separately.
    //    ed_within_k returns only entries whose distance is <= max_distance
    //    (boundary included). The pruning for this bound is exact, so no
    //    candidate within the threshold is missed.
    // ============================================================
    {
        print_header("DNA / AffineGapCost(m=1, open=2, extend=1) / ed_within_k");
        std::vector<std::string> targets = {"ACGT", "ACGTACGT", "AC", "ACGTT", "TTTT"};
        PatriciaTree<DnaAlphabet> tree(targets);
        DTPatricia<DnaAlphabet, AffineGapCost> aligner(tree, AffineGapCost(1, 2, 1));

        const std::string query = "ACGT";
        const int max_distance = 5;  // return only entries within this edit distance
        print_results(targets, query, aligner.ed_within_k(query, max_distance));
    }

    // ============================================================
    // 4) R/Y alphabet + UnitCost
    //    Collapse A,G -> R (purine) / C,T,U -> Y (pyrimidine) before
    //    comparing. "GCAT", which differs as DNA, maps to the same
    //    "ACGT" (= R Y R Y) under R/Y and gets score=0. Changing the
    //    alphabet changes the very meaning of the distance.
    // ============================================================
    {
        print_header("R/Y / UnitCost / ed_to_all  (A,G->R / C,T->Y)");
        std::vector<std::string> targets = {"ACGT", "AGGT", "GCAT", "ACGTACGT"};
        PatriciaTree<RyAlphabet> tree(targets);
        DTPatricia<RyAlphabet, UnitCost> aligner(tree);

        const std::string query = "ACGT";  // -> R Y R Y
        print_results(targets, query, aligner.ed_to_all(query));
    }

    // ============================================================
    // 5) Protein alphabet + UnitCost
    //    Handles the 20 amino acids. ed_kth_smallest collects up to top_k
    //    entries in ascending order of distance (nearest top-k). Note that
    //    when there are ties it may return more than top_k entries (see the
    //    explanation of ties in section 7).
    // ============================================================
    {
        print_header("Protein / UnitCost / ed_kth_smallest");
        std::vector<std::string> targets = {"MKVLAA", "MKVLAG", "MKILAA", "ACDEFG", "MKVL"};
        PatriciaTree<ProteinAlphabet> tree(targets);
        DTPatricia<ProteinAlphabet, UnitCost> aligner(tree);

        const std::string query = "MKVLAA";
        const std::size_t top_k = 3;  // the nearest k entries
        print_results(targets, query, aligner.ed_kth_smallest(query, top_k));
    }

    // ============================================================
    // 6) Calling search_kernel directly: search with a timeout
    //    ed_to_all and friends are thin wrappers that pass a stop predicate
    //    to search_kernel. The predicate signature is
    //        bool(int current_score, const std::vector<AlignmentResult>& results)
    //    and it is called once per edit-distance level, *after* all matches
    //    at that distance have been pushed into `results`. Returning true
    //    aborts the search. current_score is the edit distance just
    //    processed (monotonically increasing).
    //    Here we abort based on elapsed time. Because the predicate is only
    //    called per distance level, the timeout granularity is "per one
    //    edit-distance worth of work".
    // ============================================================
    {
        print_header("DNA / UnitCost / search_kernel : timeout");
        std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG", "TGCA", "GGGG"};
        PatriciaTree<DnaAlphabet> tree(targets);
        DTPatricia<DnaAlphabet, UnitCost> aligner(tree);

        const std::string query = "ACGT";
        const auto time_budget = std::chrono::milliseconds(50);
        const auto deadline = std::chrono::steady_clock::now() + time_budget;

        // Abort once the deadline passes. The dictionary is small, so in
        // practice everything is found before the deadline and the search
        // ends naturally; the goal here is just to show the mechanism.
        auto results = aligner.search_kernel(
            query, [deadline](int /*current_score*/, const std::vector<AlignmentResult> & /*r*/) {
                return std::chrono::steady_clock::now() >= deadline;
            });
        print_results(targets, query, results);
    }

    // ============================================================
    // 7) Calling search_kernel directly: top-k search with a threshold
    //    Stop at "distance <= max_dist" and "at most k entries". Passing
    //    max_dist as the 3rd argument upper_bound enables pruning and speeds
    //    things up. The pruning is exact with respect to the given bound
    //    (it never drops a candidate that could be within the bound, and the
    //    boundary value = the bound is included). upper_bound = -1 (default)
    //    disables pruning.
    //
    //    [Handling ties]
    //    The search pushes *all* matches at a given edit distance before
    //    evaluating the predicate. So r.size() >= k becomes true not at
    //    "exactly k entries" but at "having just finished pushing the
    //    edit-distance level that reached k". In other words, if several
    //    entries share the same distance as the k-th one, they are all
    //    included and the return value can exceed k. Note that truncating
    //    strictly to k would make the choice among ties arbitrary.
    // ============================================================
    {
        print_header("DNA / UnitCost / search_kernel : top-k with threshold");
        std::vector<std::string> targets = {"ACGT", "ACGA", "AAGT", "ACG",
                                            "TGCA", "GGGG", "ACGG", "ACGC"};
        PatriciaTree<DnaAlphabet> tree(targets);
        DTPatricia<DnaAlphabet, UnitCost> aligner(tree);

        const std::string query = "ACGT";
        const std::size_t k = 3;
        const int max_dist = 2;

        auto results = aligner.search_kernel(
            query,
            [k, max_dist](int current_score, const std::vector<AlignmentResult> &r) {
                return r.size() >= k || current_score >= max_dist;
            },
            max_dist);  // upper_bound = max_dist for pruning
        print_results(targets, query, results);
    }

    return 0;
}
