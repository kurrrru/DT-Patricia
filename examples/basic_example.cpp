// 後で整備された利用例に置き換える
#include <iostream>
#include <dt_patricia/dt_patricia.hpp>

int main() {
    std::vector<std::string> targets = {"ACGT", "ACG", "AC", "A", "CGT", "CG", "C", "GT", "G", "T"};
    dt_patricia::PatriciaTree tree(targets);
    dt_patricia::DTPatricia<dt_patricia::DnaAlphabet, dt_patricia::UnitCost> dtp_aligner(tree);

    std::string query = "ACG";
    auto results = dtp_aligner.ed_to_all(query);
    for (const auto& result : results) {
        std::cout << "String ID: " << result.string_id
            << ", String: " << targets[result.string_id]
            << ", Edit Distance: " << result.score << std::endl;
    }
    return 0;
}
