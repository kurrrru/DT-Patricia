# brute_test を「assert 有効（通常版）」と「assert 無効（NDEBUG 版）」の2系統でビルドする。
# 実行は tests/ をカレントにする（test.cpp が相対パス testcase/ を参照するため）。

option(DT_PATRICIA_SANITIZE
       "Build tests with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)

set(_dtp_sanitize)
if(DT_PATRICIA_SANITIZE)
  set(_dtp_sanitize
      -fsanitize=address,undefined
      -fno-sanitize-recover=all   # UBSan 検出時に即異常終了させ、ctest を確実に落とす
      -fno-omit-frame-pointer     # スタックトレースを読みやすく
      -g)
endif()

set(_dtp_warn -Wall -Wextra -Wshadow -Wold-style-cast -Wnull-dereference -Wdouble-promotion -Wformat=2 -Werror -pedantic)

# 通常版（assert 有効）
add_executable(brute_test test.cpp)
target_link_libraries(brute_test PRIVATE dt_patricia::dt_patricia)
target_compile_options(brute_test PRIVATE -O3 -UNDEBUG ${_dtp_warn} ${_dtp_sanitize})
target_link_options(brute_test    PRIVATE ${_dtp_sanitize})

# NDEBUG 版（assert 無効）
add_executable(brute_test_ndebug test.cpp)
target_link_libraries(brute_test_ndebug PRIVATE dt_patricia::dt_patricia)
target_compile_options(brute_test_ndebug PRIVATE -O2 -DNDEBUG ${_dtp_warn} ${_dtp_sanitize})
target_link_options(brute_test_ndebug    PRIVATE ${_dtp_sanitize})

# empty test
add_executable(empty_test test_empty.cpp)
target_link_libraries(empty_test PRIVATE dt_patricia::dt_patricia)
target_compile_options(empty_test PRIVATE -O2 -DNDEBUG ${_dtp_warn} ${_dtp_sanitize})
target_link_options(empty_test    PRIVATE ${_dtp_sanitize})

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
