#pragma once
#include <span>
#include <string_view>

#include <dt_patricia/patricia_tree.hpp>

namespace dt_patricia::debug {

template <AlphabetPolicy Alphabet>
std::span<const uint32_t> exact_match(const PatriciaTree<Alphabet> &tree, std::string_view query) {
    uint32_t current_node = tree.root_id();
    size_t query_pos = 0;

    while (true) {
        std::string_view label = tree.get_label(current_node);
        size_t label_len = label.length();

        size_t label_pos =
            std::mismatch(label.begin(), label.end(), query.begin() + query_pos, query.end())
                .first -
            label.begin();
        if (label_pos < label_len) {
            return {};  // 不一致
        }
        query_pos += label_len;

        if (query_pos == query.length()) {
            uint32_t leaf_node = tree.transition(current_node, PatriciaTree<Alphabet>::CODE_TERM);
            if (leaf_node != 0) {
                return tree.get_string_id(leaf_node);
            } else {
                return {};
            }
        }

        char next_char = query[query_pos];
        uint32_t next_node = tree.transition(current_node, next_char);
        if (next_node == 0) {
            return {};
        }
        current_node = next_node;
    }
}

}  // namespace dt_patricia::debug
