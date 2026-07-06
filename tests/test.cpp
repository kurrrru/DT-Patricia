// ============================================================
// test.cpp  —  Comparison test between DT-Patricia and BruteForceChecker
// ============================================================
//
// [Usage]
//   # Test all .txt files in the testcase/ directory
//   ./brute_test
//
//   # Test specific testcase files (multiple allowed)
//   ./brute_test testcase/01_unit_basic.txt testcase/04_linear_basic.txt
//
// [Testcase file format]
//   - Lines beginning with # are treated as comments and ignored
//   - Empty lines are ignored
//   - Section headers: TARGETS / QUERIES / TESTS
//   - Alphabet selection: ALPHABET <type>   (defaults to DnaAlphabet if omitted)
//     e.g.)
//       ALPHABET DnaAlphabet      — A/C/G/T
//       ALPHABET RyAlphabet       — purine(R) / pyrimidine(Y)
//       ALPHABET ProteinAlphabet  — 20 amino acids
//   - Cost selection: COST <type> [params...]
//     e.g.)
//       COST UnitCost
//       COST LinearGapCost <mismatch> <gap>
//       COST AffineGapCost <mismatch> <gap_open> <gap_extend>
//         (affine gap cost formula: gap_open + gap_extend * L)
//
//   TARGETS section: list the target strings, one per line
//   QUERIES section: list the query strings, one per line
//     (the usable characters differ per ALPHABET; DnaAlphabet allows only ATCG)
//   TESTS section: list the operations to run, one per line
//     ed_to_all              — compute the edit distance to every target
//     ed_within_k <k>        — list targets whose edit distance is <= k
//     ed_kth_smallest <k>    — list all targets with a score <= the k-th
//                              smallest score (ties included)
//
// [Testcase file example]
//   # Basic test for UnitCost
//   COST UnitCost
//   TARGETS
//   ACGT
//   ATCG
//   QUERIES
//   ACG
//   TESTS
//   ed_within_k 2
//   ed_kth_smallest 3
//   ed_to_all
//
// ============================================================

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dt_patricia/dt_patricia.hpp"

// ============================================================
// ANSI color codes
// ============================================================
static const char *COLOR_GREEN = "\033[32m";
static const char *COLOR_RED = "\033[31m";
static const char *COLOR_RESET = "\033[0m";

// ============================================================
// BruteForceChecker
// Computes the edit distance with a naive DP. Always returns the correct answer.
// ============================================================
template <typename Alphabet = dt_patricia::DnaAlphabet, typename CostType = dt_patricia::UnitCost>
class BruteForceChecker {
 public:
    BruteForceChecker(const std::vector<std::string> &targets, CostType cost = CostType())
        : _targets(targets), _cost(cost) {}

    // Return every target ordered by (score ASC, string_id ASC)
    std::vector<dt_patricia::AlignmentResult> ed_to_all(const std::string &query) const {
        std::vector<dt_patricia::AlignmentResult> results;
        results.reserve(_targets.size());
        for (uint32_t i = 0; i < static_cast<uint32_t>(_targets.size()); ++i) {
            results.push_back({i, compute_ed(_targets[i], query)});
        }
        sort_results(results);
        return results;
    }

    // Return targets with score <= k, ordered by (score ASC, string_id ASC)
    std::vector<dt_patricia::AlignmentResult> ed_within_k(const std::string &query, int k) const {
        std::vector<dt_patricia::AlignmentResult> results;
        for (uint32_t i = 0; i < static_cast<uint32_t>(_targets.size()); ++i) {
            uint32_t score = compute_ed(_targets[i], query);
            if (static_cast<int>(score) <= k) {
                results.push_back({i, score});
            }
        }
        sort_results(results);
        return results;
    }

    // Return all targets with a score <= the k-th smallest score, ordered by
    // (score ASC, string_id ASC) (ties at the boundary score are included)
    std::vector<dt_patricia::AlignmentResult> ed_kth_smallest(const std::string &query,
                                                              size_t k) const {
        if (k == 0) {
            return {};
        }
        auto all = ed_to_all(query);
        if (k >= all.size()) {
            return all;
        }
        uint32_t threshold = all[k - 1].score;
        auto it = std::upper_bound(
            all.begin(), all.end(), threshold,
            [](uint32_t val, const dt_patricia::AlignmentResult &r) { return val < r.score; });
        all.erase(it, all.end());
        return all;
    }

 private:
    static constexpr uint32_t INF = 1'000'000'000u;

    // Character matching is performed through the alphabet's code conversion.
    // This way, even when several characters collapse onto the same code (as in
    // RyAlphabet), matching is identical to DTPatricia, making this a correct oracle.
    static constexpr std::array<uint8_t, 256> CHAR_TO_CODE = Alphabet::make_char_to_code();
    static uint8_t code_of(char c) { return CHAR_TO_CODE[static_cast<uint8_t>(c)]; }

    static void sort_results(std::vector<dt_patricia::AlignmentResult> &v) {
        std::sort(v.begin(), v.end(),
                  [](const dt_patricia::AlignmentResult &a, const dt_patricia::AlignmentResult &b) {
                      return a.score != b.score ? a.score < b.score : a.string_id < b.string_id;
                  });
    }

    uint32_t compute_ed(const std::string &s, const std::string &t) const {
        if constexpr (CostType::is_linear) {
            return compute_ed_linear(s, t);
        } else {
            return compute_ed_affine(s, t);
        }
    }

    // Standard edit-distance DP (UnitCost / LinearGapCost)
    uint32_t compute_ed_linear(const std::string &s, const std::string &t) const {
        const size_t n = s.size(), m = t.size();
        std::vector<std::vector<uint32_t>> dp(n + 1, std::vector<uint32_t>(m + 1, INF));
        dp[0][0] = 0;
        for (size_t i = 1; i <= n; ++i) {
            dp[i][0] = _cost.gap * static_cast<uint32_t>(i);
        }
        for (size_t j = 1; j <= m; ++j) {
            dp[0][j] = _cost.gap * static_cast<uint32_t>(j);
        }
        for (size_t i = 1; i <= n; ++i) {
            for (size_t j = 1; j <= m; ++j) {
                uint32_t sub = (code_of(s[i - 1]) == code_of(t[j - 1])) ? 0u : _cost.mismatch;
                dp[i][j] = std::min(
                    {dp[i - 1][j - 1] + sub, dp[i - 1][j] + _cost.gap, dp[i][j - 1] + _cost.gap});
            }
        }
        return dp[n][m];
    }

    // Gotoh's algorithm (AffineGapCost)
    // cost(L) = gap_open + gap_extend * L
    uint32_t compute_ed_affine(const std::string &s, const std::string &t) const {
        const size_t n = s.size(), m = t.size();
        // H[i][j]: optimal cost for s[0..i-1] and t[0..j-1]
        // E[i][j]: optimal cost ending with a gap on the t side
        // F[i][j]: optimal cost ending with a gap on the s side
        std::vector<std::vector<uint32_t>> H(n + 1, std::vector<uint32_t>(m + 1, INF));
        std::vector<std::vector<uint32_t>> E(n + 1, std::vector<uint32_t>(m + 1, INF));
        std::vector<std::vector<uint32_t>> F(n + 1, std::vector<uint32_t>(m + 1, INF));

        H[0][0] = 0;
        for (size_t i = 1; i <= n; ++i) {
            H[i][0] = _cost.gap_open + _cost.gap_extend * static_cast<uint32_t>(i);
            F[i][0] = H[i][0];
        }
        for (size_t j = 1; j <= m; ++j) {
            H[0][j] = _cost.gap_open + _cost.gap_extend * static_cast<uint32_t>(j);
            E[0][j] = H[0][j];
        }
        for (size_t i = 1; i <= n; ++i) {
            for (size_t j = 1; j <= m; ++j) {
                E[i][j] = std::min(H[i][j - 1] + _cost.gap_open + _cost.gap_extend,
                                   E[i][j - 1] + _cost.gap_extend);
                F[i][j] = std::min(H[i - 1][j] + _cost.gap_open + _cost.gap_extend,
                                   F[i - 1][j] + _cost.gap_extend);
                uint32_t sub = (code_of(s[i - 1]) == code_of(t[j - 1])) ? 0u : _cost.mismatch;
                H[i][j] = std::min({H[i - 1][j - 1] + sub, E[i][j], F[i][j]});
            }
        }
        return H[n][m];
    }

    std::vector<std::string> _targets;
    CostType _cost;
};

// ============================================================
// Testcase struct and parser
// ============================================================

struct TestOp {
    std::string name;
    int k = 0;
};

struct TestCase {
    std::string filename;
    std::string alphabet = "DnaAlphabet";  // default when omitted
    std::string cost_type;
    uint32_t cost_p1 = 1, cost_p2 = 1, cost_p3 = 1;
    std::vector<std::string> targets;
    std::vector<std::string> queries;
    std::vector<TestOp> ops;
};

// Set of characters allowed for each alphabet (uppercase only).
// Returns an empty string for an unknown alphabet name.
static std::string allowed_chars(const std::string &alphabet) {
    if (alphabet == "DnaAlphabet") {
        return "ACGT";
    }
    if (alphabet == "RyAlphabet") {
        return "ACGTURY";
    }
    if (alphabet == "ProteinAlphabet") {
        return "ACDEFGHIKLMNPQRSTVWY";
    }
    return "";
}

static bool is_valid_seq(const std::string &s, const std::string &allowed) {
    if (s.empty()) {
        return false;
    }
    for (char c : s) {
        if (allowed.find(c) == std::string::npos) {
            return false;
        }
    }
    return true;
}

static TestCase parse_test_file(const std::string &path) {
    TestCase tc;
    tc.filename = path;
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open: " + path);
    }

    std::string line, section;
    while (std::getline(file, line)) {
        while (!line.empty() && (line.back() == ' ' || line.back() == '\r')) {
            line.pop_back();
        }
        if (line.empty() || line[0] == '#') {
            continue;
        }

        if (line == "TARGETS" || line == "QUERIES" || line == "TESTS") {
            section = line;
            continue;
        }
        if (line.rfind("ALPHABET ", 0) == 0) {
            std::istringstream iss(line.substr(9));
            iss >> tc.alphabet;
            if (allowed_chars(tc.alphabet).empty()) {
                throw std::runtime_error("Unknown alphabet: \"" + tc.alphabet + "\"");
            }
            continue;
        }
        if (line.rfind("COST ", 0) == 0) {
            std::istringstream iss(line.substr(5));
            iss >> tc.cost_type;
            if (tc.cost_type == "LinearGapCost") {
                iss >> tc.cost_p1 >> tc.cost_p2;
            }
            if (tc.cost_type == "AffineGapCost") {
                iss >> tc.cost_p1 >> tc.cost_p2 >> tc.cost_p3;
            }
            continue;
        }

        const std::string allowed = allowed_chars(tc.alphabet);
        if (section == "TARGETS") {
            if (!is_valid_seq(line, allowed)) {
                throw std::runtime_error("Invalid target for " + tc.alphabet + ": \"" + line +
                                         "\"");
            }
            tc.targets.push_back(line);
        } else if (section == "QUERIES") {
            if (!is_valid_seq(line, allowed)) {
                throw std::runtime_error("Invalid query for " + tc.alphabet + ": \"" + line + "\"");
            }
            tc.queries.push_back(line);
        } else if (section == "TESTS") {
            TestOp op;
            std::istringstream iss(line);
            iss >> op.name;
            if (op.name == "ed_within_k" || op.name == "ed_kth_smallest") {
                iss >> op.k;
            }
            tc.ops.push_back(op);
        }
    }

    if (tc.targets.empty()) {
        throw std::runtime_error("No targets in: " + path);
    }
    if (tc.queries.empty()) {
        throw std::runtime_error("No queries in: " + path);
    }
    if (tc.ops.empty()) {
        throw std::runtime_error("No tests in: " + path);
    }
    return tc;
}

// ============================================================
// Result comparison utilities
// ============================================================

static void sort_results(std::vector<dt_patricia::AlignmentResult> &v) {
    std::sort(v.begin(), v.end(),
              [](const dt_patricia::AlignmentResult &a, const dt_patricia::AlignmentResult &b) {
                  return a.score != b.score ? a.score < b.score : a.string_id < b.string_id;
              });
}

static bool results_equal(const std::vector<dt_patricia::AlignmentResult> &a,
                          const std::vector<dt_patricia::AlignmentResult> &b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].string_id != b[i].string_id || a[i].score != b[i].score) {
            return false;
        }
    }
    return true;
}

static void print_mismatch(const std::string &query, const std::string &op_name, int k,
                           const std::vector<dt_patricia::AlignmentResult> &expected,
                           const std::vector<dt_patricia::AlignmentResult> &actual,
                           const std::vector<std::string> &targets) {
    std::cout << "    MISMATCH query=\"" << query << "\" op=" << op_name;
    if (op_name != "ed_to_all") {
        std::cout << "(" << k << ")";
    }
    std::cout << "\n";

    std::cout << "    brute [" << expected.size() << "]:";
    for (const auto &r : expected) {
        std::cout << " " << targets[r.string_id] << ":" << r.score;
    }
    std::cout << "\n";

    std::cout << "    dtp   [" << actual.size() << "]:";
    for (const auto &r : actual) {
        std::cout << " " << targets[r.string_id] << ":" << r.score;
    }
    std::cout << "\n";
}

// ============================================================
// Test execution
// ============================================================

template <typename Alphabet, typename CostType>
static bool run_with_cost(const TestCase &tc, CostType cost) {
    dt_patricia::PatriciaTree<Alphabet> tree(tc.targets);
    dt_patricia::DTPatricia<Alphabet, CostType> dtp_aligner(tree, cost);
    BruteForceChecker<Alphabet, CostType> brute(tc.targets, cost);

    bool all_passed = true;
    for (const auto &query : tc.queries) {
        for (const auto &op : tc.ops) {
            std::vector<dt_patricia::AlignmentResult> expected, actual;

            if (op.name == "ed_to_all") {
                expected = brute.ed_to_all(query);
                actual = dtp_aligner.ed_to_all(query);
            } else if (op.name == "ed_within_k") {
                expected = brute.ed_within_k(query, op.k);
                actual = dtp_aligner.ed_within_k(query, op.k);
            } else if (op.name == "ed_kth_smallest") {
                expected = brute.ed_kth_smallest(query, static_cast<size_t>(op.k));
                actual = dtp_aligner.ed_kth_smallest(query, static_cast<size_t>(op.k));
            } else {
                std::cout << "    Unknown op: " << op.name << "\n";
                all_passed = false;
                continue;
            }

            sort_results(expected);
            sort_results(actual);

            if (!results_equal(expected, actual)) {
                print_mismatch(query, op.name, op.k, expected, actual, tc.targets);
                all_passed = false;
            }
        }
    }
    return all_passed;
}

template <typename Alphabet>
static bool dispatch_cost(const TestCase &tc) {
    if (tc.cost_type == "UnitCost") {
        return run_with_cost<Alphabet>(tc, dt_patricia::UnitCost{});
    }
    if (tc.cost_type == "LinearGapCost") {
        return run_with_cost<Alphabet>(tc, dt_patricia::LinearGapCost(tc.cost_p1, tc.cost_p2));
    }
    if (tc.cost_type == "AffineGapCost") {
        return run_with_cost<Alphabet>(
            tc, dt_patricia::AffineGapCost(tc.cost_p1, tc.cost_p2, tc.cost_p3));
    }

    std::cout << "  Unknown cost type: " << tc.cost_type << "\n";
    return false;
}

static bool run_test_case(const TestCase &tc) {
    if (tc.alphabet == "DnaAlphabet") {
        return dispatch_cost<dt_patricia::DnaAlphabet>(tc);
    }
    if (tc.alphabet == "RyAlphabet") {
        return dispatch_cost<dt_patricia::RyAlphabet>(tc);
    }
    if (tc.alphabet == "ProteinAlphabet") {
        return dispatch_cost<dt_patricia::ProteinAlphabet>(tc);
    }

    std::cout << "  Unknown alphabet: " << tc.alphabet << "\n";
    return false;
}

static bool add_files(const std::string &path, std::vector<std::string> &files) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!fs::exists(path, ec)) {
        if (ec) {
            std::cerr << "Error checking existence of file: " << path << " (" << ec.message()
                      << ")\n";
        } else {
            std::cerr << "Error: file not found: " << path << "\n";
        }
        return false;
    }
    if (fs::is_directory(path, ec)) {
        for (const auto &entry : fs::directory_iterator(path, ec)) {
            if (entry.path().extension() == ".txt") {
                files.push_back(entry.path().string());
            }
        }
    } else if (fs::is_regular_file(path, ec)) {
        files.push_back(path);
    } else {
        std::cerr << "Error: not a regular file or directory: " << path << "\n";
        return false;
    }
    return true;
}

// ============================================================
// main
// ============================================================

int main(int argc, char *argv[]) {
    namespace fs = std::filesystem;

    // Collect the list of files to run
    std::vector<std::string> files;
    bool has_error = false;

    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            has_error |= !add_files(argv[i], files);
        }
    } else {
        has_error |= !add_files("testcase", files);
    }
    if (files.empty() && !(argc == 2 && std::string(argv[1]) == "testcase_random")) {
        std::cerr << "No .txt files to run.\n";
        return 1;
    }
    std::sort(files.begin(), files.end());

    int passed = 0, failed = 0;
    for (const auto &f : files) {
        bool ok = false;
        std::string error_msg;
        try {
            TestCase tc = parse_test_file(f);
            ok = run_test_case(tc);
        } catch (const std::exception &e) {
            error_msg = e.what();
        }

        if (!error_msg.empty()) {
            std::cout << COLOR_RED << "[ERROR]" << COLOR_RESET << " " << f << "\n"
                      << "    " << error_msg << "\n";
            ++failed;
        } else if (ok) {
            std::cout << COLOR_GREEN << "[PASSED]" << COLOR_RESET << " " << f << "\n";
            ++passed;
        } else {
            std::cout << COLOR_RED << "[FAILED]" << COLOR_RESET << " " << f << "\n";
            ++failed;
        }
    }

    std::cout << "\n" << passed << "/" << (passed + failed) << " passed\n";
    return (has_error || failed > 0) ? 1 : 0;
}
