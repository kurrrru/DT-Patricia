ドキュメント整備時に書き直すので、今はメモだけ

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


# テスト
```sh
# === テスト（正しさ）===
cmake -S . -B build/tests \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DDTPB_BUILD_TESTS=ON \
  -DDTPB_BUILD_BENCHMARK=OFF
cmake --build build/tests
ctest --test-dir build/tests --output-on-failure
```


```sh
# === テスト（正しさ）===
cmake --preset tests          # configure → build/tests/ に生成
cmake --build --preset tests  # ビルド
ctest --preset tests          # 実行（--output-on-failure は preset に内蔵）
```

```sh
# 生成（既定: シード42, 60件, tests/testcase_random へ）
python3 tests/scripts/gen_random_tests.py
python3 tests/scripts/gen_random_tests.py --seed 42 --count 300   # 変更も可
```
