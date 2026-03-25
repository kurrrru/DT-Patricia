#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <cstring>

// PatriciaWFA
#include "PatriciaTree.hpp"
#include "PatriciaWFA.hpp"
#include "CostType.hpp"

// ---------------------------------------------------------
// Utility: Timer
// ---------------------------------------------------------
class ScopedTimer {
    using Clock = std::chrono::high_resolution_clock;
    
    std::string _name;
    Clock::time_point _start;
    long long _duration_sum;

    static constexpr Clock::time_point kInvalidTime = Clock::time_point::min();

public:
    ScopedTimer(std::string name) 
        : _name(name), _start(kInvalidTime), _duration_sum(0) {}

    void start() {
        auto now_time = Clock::now();

        if (_start != kInvalidTime) {
            std::cerr << "[WARNING] " << _name << ": start() called twice without end()!" << std::endl;
        }

        _start = now_time;
    }

    void end() {
        auto end_time = Clock::now();

        if (_start == kInvalidTime) {
            std::cerr << "[WARNING] " << _name << ": end() called without start()!" << std::endl;
            return; 
        }

        auto ms = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - _start).count();
        _duration_sum += ms;

        _start = kInvalidTime;
    }

    ~ScopedTimer() {
        if (_start != kInvalidTime) {
            std::cerr << "[WARNING] " << _name << ": Destructor called while timer is still running! (forgot end()?)" << std::endl;
        }

        double us = _duration_sum / 1000.0;
        std::cout << "[TIME] " << std::left << std::setw(30) << _name << ": " << std::fixed << std::setprecision(4) <<
        us << " us" << std::endl;
    }
};

// ---------------------------------------------------------
// Utility: Data Loader
// ---------------------------------------------------------
std::vector<std::string> load_lines_fast(const std::string& filepath) {
    std::ifstream f(filepath, std::ios::binary);
    if (!f.is_open()) {
        std::cerr << "Error: Could not open " << filepath << std::endl;
        exit(1);
    }
    
    // Read entire file
    f.seekg(0, std::ios::end);
    size_t file_size = f.tellg();
    f.seekg(0, std::ios::beg);
    
    std::string buffer(file_size, '\0');
    f.read(&buffer[0], file_size);
    
    // Split by newline
    std::vector<std::string> lines;
    lines.reserve(file_size / 100);
    
    size_t start = 0;
    size_t line_num = 1;

    for (size_t i = 0; i < file_size; ++i) {
        if (buffer[i] == '\n') {
            if (i > start) {
                std::string line = buffer.substr(start, i - start);

                bool invalid_char_found = false;
                for (size_t k = 0; k < line.size(); ++k) {
                    char c = line[k];
                    if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                        std::cerr << "[WARN] Invalid char '" << c << "' (code: " << (int)c 
                                  << ") found in " << filepath << " at line " << line_num 
                                  << ", pos " << k+1 << std::endl;
                        invalid_char_found = true;
                    }
                }

                if (!invalid_char_found) {
                    lines.emplace_back(std::move(line));
                }
            }
            start = i + 1;
            line_num++;
        }
    }
    
    // Last line
    if (start < file_size) {
        std::string line = buffer.substr(start);
        
        bool invalid_char_found = false;
        for (size_t k = 0; k < line.size(); ++k) {
            char c = line[k];
            if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                std::cerr << "[WARN] Invalid char '" << c << "' (code: " << (int)c 
                          << ") found in " << filepath << " at last line" << std::endl;
                invalid_char_found = true;
            }
        }
        
        if (!invalid_char_found) {
            lines.emplace_back(std::move(line));
        }
    }
    
    return lines;
}

// ---------------------------------------------------------
// Main Benchmark Routine
// ---------------------------------------------------------
void print_usage(const char* prog_name) {
    std::cerr << "Usage: " << prog_name << " <db_file> <query_file> <mode>\n";
    std::cerr << "Modes:\n";
    std::cerr << "  range    : Run ed_within_k only (Patricia)\n";
    std::cerr << "  standard : Run ed_within_k + ed_kth_smallest (Patricia)\n";
    std::cerr << "  full     : Run ALL (ed_within_k, ed_kth_smallest, ed_to_all)\n";
}

int main(int argc, char** argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string db_path = argv[1];
    std::string query_path = argv[2];
    std::string mode = argv[3]; // range, standard, full

    std::cout << "=======================================================\n";
    std::cout << "BENCHMARK (Affine Gap): " << db_path << "\n";
    std::cout << "MODE: " << mode << "\n";
    std::cout << "=======================================================\n";

    // Load data
    auto dictionary = load_lines_fast(db_path);
    auto queries = load_lines_fast(query_path);

    std::cout << "Dictionary: " << dictionary.size() << "\n";
    std::cout << "Queries: " << queries.size() << "\n";

    // Parameter sets: (mismatch, gap_open, gap_extend)
    std::vector<std::tuple<int, int, int>> param_sets = {
        {4, 6, 2},  // BWA-MEM default
        {6, 5, 3}   // Bowtie2 default
    };

    for (const auto& [mismatch, gap_open, gap_extend] : param_sets) {
        std::cout << "\n=======================================================\n";
        std::cout << "Parameters: (mismatch, gap_open, gap_extend) = (" 
                  << mismatch << ", " << gap_open << ", " << gap_extend << ")\n";
        std::cout << "=======================================================\n";

        // Build Trie
        std::cout << "\n[Building PatriciaTrie]..." << std::endl;
        PatriciaTree pat_trie;
        {
            ScopedTimer t("PatriciaTrie Build");
            t.start();
            pat_trie = PatriciaTree(dictionary);
            t.end();
        }

        // Prepare aligners
        AffineGapCost affine_cost(mismatch, gap_open, gap_extend);
        PatriciaWFA<AffineGapCost> pat_aligner(pat_trie, affine_cost);

        // -------------------------------------------------
        // 1. Range Query (ed_within_k) - Patricia only
        // -------------------------------------------------
        std::cout << "\n[Bench: ed_within_k (Range Query)]" << std::endl;
        
        // Base K values (for edit distance)
        std::vector<int> K_VALS_BASE = {1, 3, 5, 10, 20, 50, 100};
        
        // Scale K values for affine gap based on mismatch
        std::vector<int> K_VALS_SCALED;
        for (int k : K_VALS_BASE) {
            K_VALS_SCALED.push_back(k * mismatch);
        }
        
        for (size_t i = 0; i < K_VALS_SCALED.size(); ++i) {
            int k = K_VALS_SCALED[i];
            std::cout << "--- k = " << k << " (base k = " << K_VALS_BASE[i] << ") ---" << std::endl;
            
            // Patricia
            {
                ScopedTimer t("PatriciaWFA");
                volatile size_t checksum = 0;
                for(const auto& q : queries) {
                    t.start();
                    auto results = pat_aligner.ed_within_k(q, k);
                    t.end();
                    checksum += results.size();
                }
                std::cout << " Checksum: " << checksum << std::endl;
            }

            // // KSW2
            // {
            //     ScopedTimer t("KSW2 (Global Affine)");
            //     Ksw2AffineBench ksw2_bench(dictionary, mismatch, gap_open, gap_extend);
            //     volatile long long checksum = 0;
            //     for(const auto& q : queries) {
            //         checksum += ksw2_bench.ed_within_k(q, k, t);
            //     }
            //     std::cout << " Checksum: " << checksum << std::endl;
            // }

            // // WFA2-lib
            // {
            //     ScopedTimer t("WFA2-lib (Affine)");
            //     Wfa2AffineBench wfa2_bench(dictionary, mismatch, gap_open, gap_extend);
            //     volatile long long checksum = 0;
            //     for(const auto& q : queries) {
            //         checksum += wfa2_bench.ed_within_k(q, k, t);
            //     }
            //     std::cout << " Checksum: " << checksum << std::endl;
            // }
        }

        if (mode == "range") continue; // Next parameter set

        // -------------------------------------------------
        // 2. Top-k Query (ed_kth_smallest) - Patricia only
        // -------------------------------------------------
        std::cout << "\n[Bench: ed_kth_smallest (Top-k Query)]" << std::endl;
        std::cout << "(Skipped in this benchmark to focus on ed_within_k and ed_to_all)" << std::endl;
        // std::vector<int> TOP_K_VALS = {1, 5, 10, 20, 50, 100};

        // for (int k : TOP_K_VALS) {
        //     std::cout << "--- Top-" << k << " ---" << std::endl;
            
        //     // Patricia
        //     {
        //         ScopedTimer t("PatriciaWFA");
        //         volatile size_t checksum = 0;
        //         for(const auto& q : queries) {
        //             t.start();
        //             auto results = pat_aligner.ed_kth_smallest(q, k);
        //             t.end();
        //             checksum += results.size();
        //         }
        //         std::cout << " Checksum: " << checksum << std::endl;
        //     }
        // }

        if (mode == "standard") continue; // Next parameter set

        // -------------------------------------------------
        // 3. Full Scan (ed_to_all) - All algorithms
        // -------------------------------------------------
        std::cout << "\n[Bench: ed_to_all (Full Scan)]" << std::endl;
        
        // Patricia
        // {
        //     ScopedTimer t("PatriciaWFA");
        //     volatile size_t checksum = 0;
        //     for(const auto& q : queries) {
        //         t.start();
        //         auto results = pat_aligner.ed_to_all(q);
        //         t.end();
        //         checksum += results.size(); 
        //     }
        //     std::cout << " Checksum: " << checksum << std::endl;
        // }
        
        // WFA2-lib
        // {
        //     ScopedTimer t("WFA2-lib (Affine)");
        //     Wfa2AffineBench wfa2_bench(dictionary, mismatch, gap_open, gap_extend);
        //     volatile long long checksum = 0;
        //     for(const auto& q : queries) {
        //         checksum += wfa2_bench.bench_to_all(q, t);
        //     }
        //     std::cout << " Checksum: " << checksum << std::endl;
        // }
        
        // // KSW2
        // if (queries.size() * dictionary.size() * queries[0].size() * dictionary[0].size() <= 1e9) {
        //     ScopedTimer t("KSW2 (Global Affine)");
        //     Ksw2AffineBench ksw2_bench(dictionary, mismatch, gap_open, gap_extend);
        //     volatile long long checksum = 0;
        //     for(const auto& q : queries) {
        //         checksum += ksw2_bench.bench_to_all(q, t);
        //     }
        //     std::cout << " Checksum: " << checksum << std::endl;
        // } else {
        //     std::cout << "KSW2 skipped for ed_to_all due to excessive computation size." << std::endl;
        // }
    }
    
    return 0;
}
