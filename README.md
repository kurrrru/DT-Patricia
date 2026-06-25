ドキュメント整備時に書き直す

exampleの実行方法
```bash
# ① configure: -S=ソース場所(.) -B=ビルド場所(build)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# ② build
cmake --build build

# ③ 実行
./build/examples/basic_example
```

もしくは

```bash
cmake --preset debug
cmake --build --preset debug
./build/debug/examples/basic_example
```
