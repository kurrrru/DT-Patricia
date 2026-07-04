# brute_test を「assert 有効（通常版）」と「assert 無効（NDEBUG 版）」の2系統でビルドする。
# 実行は tests/ をカレントにする（test.cpp が相対パス testcase/ を参照するため）。

set(_dtp_warn -Wall -Wextra -Wshadow -Wold-style-cast -Wnull-dereference -Wdouble-promotion -Wformat=2 -Werror -pedantic)

# 通常版（assert 有効）
add_executable(brute_test test.cpp)
target_link_libraries(brute_test PRIVATE dt_patricia::dt_patricia)
# 末尾の -UNDEBUG で、Release ビルド時に注入される -DNDEBUG を打ち消し assert を必ず有効化
target_compile_options(brute_test PRIVATE -O3 -UNDEBUG ${_dtp_warn})

# NDEBUG 版（assert 無効）
add_executable(brute_test_ndebug test.cpp)
target_link_libraries(brute_test_ndebug PRIVATE dt_patricia::dt_patricia)
target_compile_options(brute_test_ndebug PRIVATE -O2 -DNDEBUG ${_dtp_warn})

# empty test
add_executable(empty_test test_empty.cpp)
target_link_libraries(empty_test PRIVATE dt_patricia::dt_patricia)
target_compile_options(empty_test PRIVATE -O2 -DNDEBUG ${_dtp_warn})

add_test(NAME brute_test
    COMMAND brute_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME brute_test_ndebug
    COMMAND brute_test_ndebug
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME empty_test
    COMMAND empty_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
# random test
add_test(NAME random_test
    COMMAND brute_test testcase_random
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
