#pragma once
#include <dt_patricia/patricia_tree.hpp>
#include <iostream>
#include <string>

namespace patricia_debug {

    namespace internal {
        inline void dump_recursive(const PatriciaTree& tree, uint32_t node_id, int depth,
                                   const std::vector<uint32_t>& counts,
                                   const std::vector<uint32_t>& min_lens,
                                   const std::vector<uint32_t>& max_lens) {
            std::string indent(depth * 2, ' ');
            
            std::cout << indent << "node" << node_id << " [label=\"" 
                      << node_id << (tree.is_terminal(node_id) ? "\\n(TERM)" : "") 
                      << "(" <<  counts[node_id] << ")" 
                      << min_lens[node_id] << "/" << max_lens[node_id]
                      << "\"];\n";

            // _base や _check を直接見ず、公開APIである transition を用いて子ノードを列挙する
            for (uint8_t code = 0; code <= 5; ++code) { 
                uint32_t next_id = tree.transition(node_id, code);
                
                if (next_id != 0 && next_id != tree.root_id()) {            
                    std::string edge_str;
                    if (code == PatriciaTree::CODE_TERM) {
                        edge_str = "$"; 
                    } else {
                        edge_str = std::string(tree.get_label(next_id));
                    }
                    std::cout << indent << "node" << node_id << " -> node" << next_id
                              << " [label=\"" << edge_str << "\"];\n";
                    
                    dump_recursive(tree, next_id, depth + 1, counts, min_lens, max_lens);
                }
            }
        }
    } // namespace internal

    inline void dump_to_dot(const PatriciaTree& tree) {
        std::cout << "digraph PatriciaTree {\n";
        
        // 再帰のたびにゲッターを呼ぶオーバーヘッドを避けるため、一度だけ取得して参照渡しする
        internal::dump_recursive(tree, tree.root_id(), 0,
                                 tree.get_subtree_counts(),
                                 tree.get_subtree_min_lengths(),
                                 tree.get_subtree_max_lengths());
                                 
        std::cout << "}\n";
    }

} // namespace patricia_debug
