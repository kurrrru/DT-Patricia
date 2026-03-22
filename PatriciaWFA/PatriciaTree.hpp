#pragma once
#include <array>
#include <vector>
#include <string>
#include <string_view>
#include <span>

class PatriciaTree {
public:
    // =========================================================
    // 1. コンストラクタ / デストラクタ (Rule of Five)
    // =========================================================
    
    // デフォルトコンストラクタ
    PatriciaTree() = default;

    // コンストラクタ（内部で build を呼ぶ）
    explicit PatriciaTree(const std::vector<std::string> &input_data);

    // デストラクタ
    ~PatriciaTree() = default;

    // コピー
    // [NOTE]辞書は大きくなるので、コピーは一旦禁止しておく
    PatriciaTree(const PatriciaTree&) = delete;
    PatriciaTree &operator=(const PatriciaTree&) = delete;

    // ムーブ
    PatriciaTree(PatriciaTree&&) noexcept = default;
    PatriciaTree &operator=(PatriciaTree&&) noexcept = default;

    // =========================================================
    // 2. 公開API
    // =========================================================

    // ツリーが空かどうかを返す
    [[nodiscard]] inline bool empty() const noexcept {
        return _base.empty();
    }

    // ルートノードのIDを取得
    // ノード0を無効値とし、ルートをノード1とする
    [[nodiscard]] inline uint32_t root_id() const noexcept {
        return 1;
    }

    // ノード node_id から 文字 ch で遷移できるか？
    // 遷移できれば次のノードID、できなければ 0 (無効値) を返す
    [[nodiscard]] inline uint32_t transition(uint32_t node_id, char ch) const noexcept {
        uint8_t code = CHAR_TO_CODE[static_cast<uint8_t>(ch)];
        return transition(node_id, code);
    }

    // ノード node_id から エンコードされた文字 code で遷移できるか？
    // エンコードは、{'\0':0, 'A':1, 'C':2, 'G':3, 'T':4} のように行う
    // 遷移できれば次のノードID、できなければ 0 (無効値) を返す
    [[nodiscard]] inline uint32_t transition(uint32_t node_id, uint8_t code) const noexcept {
        uint32_t base = _base[node_id];
        if (base == 0) {  // node_id が葉の場合
            return 0;
        }
        uint32_t next_id = base + code;
        if (next_id >= _check.size() || _check[next_id] != node_id) {
            // node_id が next_id の親でない(遷移できない)場合
            return 0;
        }
        return next_id;
    }

    // ノード node_id に入ってくるエッジのラベルを取得
    [[nodiscard]] inline std::string_view get_label(uint32_t node_id) const noexcept {
        uint32_t offset = _label_offset[node_id];
        uint32_t length = _label_len[node_id];
        return std::string_view(TEXT_POOL.data() + offset, length);
    }

    // ノード node_id が葉かどうかを返す
    // _base[node_id] が 0 なら葉とみなす
    [[nodiscard]] inline bool is_leaf(uint32_t node_id) const noexcept {
        return _base[node_id] == 0;
    }

    // ノード node_id が単語の終端かどうかを返す
    // '\0'で遷移するエッジがあれば終端とみなす
    // [NOTE]単語IDは葉ノードに紐づいており、終端ノードではない
    [[nodiscard]] inline bool is_terminal(uint32_t node_id) const noexcept {
        return transition(node_id, CODE_TERM) != 0;
    }

    // ノード node_id が葉であれば、その単語IDリストを返す
    [[nodiscard]] inline std::span<const uint32_t> get_string_id(uint32_t node_id) const noexcept {
        if (!is_leaf(node_id)) {
            return std::span<const uint32_t>();
        }
        uint32_t offset = _string_ids_offset[node_id];
        uint32_t count = _string_ids_count[node_id];
        return std::span<const uint32_t>(_string_ids.data() + offset, count);
    }

    // [TODO]: string_count()という名前の方が良い
    [[nodiscard]] inline uint32_t size() const noexcept {
        return _size;
    }

    [[nodiscard]] inline uint32_t node_count() const noexcept {
        return static_cast<uint32_t>(_base.size());
    }

    [[nodiscard]] inline uint32_t get_label_length(uint32_t node_id) const noexcept {
        return _label_len[node_id];
    }

    [[nodiscard]] inline const std::vector<uint32_t>& get_subtree_counts() const noexcept {
        return _subtree_counts;
    }

    [[nodiscard]] inline const std::vector<uint32_t>& get_subtree_max_lengths() const noexcept {
        return _subtree_max_len;
    }

    [[nodiscard]] inline const std::vector<uint32_t>& get_subtree_min_lengths() const noexcept {
        return _subtree_min_len;
    }

    [[nodiscard]] uint32_t get_parent(uint32_t node_id) const noexcept {
        if (node_id >= _check.size()) {
            return 0;
        }
        return _check[node_id];
    }

private:
    uint32_t _size;  // 登録されている文字列数

    std::vector<uint32_t> _base;  // indexはノードID、値はベース値
    std::vector<uint32_t> _check;  // indexはノードID、値は親ノードID（不正ノード、ルートノードの場合は0）
    
    // ラベル情報
    // std::string_viewを使わずに、オフセットと長さで管理する方式
    // TEXT_POOLがリサイズされてもラベル情報が壊れないようにするため
    std::vector<uint32_t> _label_offset;  // indexはノードID
    std::vector<uint32_t> _label_len;  // indexはノードID

    // 文字列の実体プール
    // _label_offsetと_label_lenで参照される
    std::string TEXT_POOL; 

    // ノードIDから対応する単語IDを取得するテーブル
    // そのノードが単語終端であるかを事前に必ず確認すること
    std::vector<uint32_t> _string_ids_offset; // index=ノードID, 値=_string_idsプール内の開始位置
    std::vector<uint32_t> _string_ids_count;  // index=ノードID, 値=そのノードに紐づくIDの個数
    std::vector<uint32_t> _string_ids;        // 全てのIDを隙間なく詰め込んだ巨大配列

    std::vector<uint32_t> _subtree_counts;   // index=ノードID, 値=そのノードを根とする部分木に含まれる単語数
    std::vector<uint32_t> _subtree_max_len;  // index=ノードID, 値=そのノードから部分木の葉までの文字列の最大長(そのノードのラベル長を含む)
    std::vector<uint32_t> _subtree_min_len;  // index=ノードID, 値=そのノードから部分木の葉までの文字列の最小長(そのノードのラベル長を含む)

    // 内部ビルド関数
    struct QueueItem {
        uint32_t node_id;
        size_t start_idx;
        size_t end_idx;
        size_t label_offset;
    };
    
    void build(const std::vector<std::string> &input_data, const std::vector<uint32_t> &sorted_indices);
    void build_node_bfs(uint32_t node_id,
                        size_t start_idx,
                        size_t end_idx,
                        size_t label_offset,
                        const std::vector<std::string> &input_data,
                        const std::vector<uint32_t> &sorted_indices,
                        std::vector<QueueItem> &queue);
    void compute_subtree_length_bounds();

 public:
    // =========================================================
    // char -> code 変換テーブル
    // =========================================================

    // SIMDの境界外読み込みアクセス違反を防ぐためのパディングサイズ
    static constexpr size_t SIMD_PADDING_SIZE = 256;

    // 定数定義
    static constexpr uint8_t CODE_TERM = 0;
    static constexpr uint8_t CODE_A    = 1;
    static constexpr uint8_t CODE_C    = 2;
    static constexpr uint8_t CODE_G    = 3;
    static constexpr uint8_t CODE_T    = 4;  // UもTとして扱う
    static constexpr uint8_t CODE_N    = 5;  // Nや未知の文字用

    static constexpr uint8_t CODE_MAX  = 5;  // 最大コード値 (Nのコード)

    // テーブル生成用ラムダ（コンパイル時に計算完了）
    static constexpr std::array<uint8_t, 256> CHAR_TO_CODE = []() consteval {
        std::array<uint8_t, 256> table{};        
        table.fill(CODE_N);

        // 1. 終端文字 (\0)
        table['\0'] = CODE_TERM;

        // 2. DNA塩基 (大文字・小文字)
        table['A'] = CODE_A; table['a'] = CODE_A;
        table['C'] = CODE_C; table['c'] = CODE_C;
        table['G'] = CODE_G; table['g'] = CODE_G;
        table['T'] = CODE_T; table['t'] = CODE_T;

        // 3. RNA塩基のUもTとして扱う
        table['U'] = CODE_T; table['u'] = CODE_T;

        // 4. その他の文字はすべてCODE_Nのまま

        return table;
    }();
};
