cmake_minimum_required(VERSION 3.25)

project(dt_patricia
    VERSION 0.1.0
    DESCRIPTION "Exact sequence-to-dictionary alignment using diagonal transition on a Patricia tree"
    LANGUAGES CXX
)

# ---- header-only library ----
add_library(dt_patricia INTERFACE)
add_library(dt_patricia::dt_patricia ALIAS dt_patricia)

target_include_directories(dt_patricia INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(dt_patricia INTERFACE cxx_std_20)

# ---- examples ----
option(DT_PATRICIA_BUILD_EXAMPLES "Build the example programs" ${PROJECT_IS_TOP_LEVEL})
if(DT_PATRICIA_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# ---- install ----
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS dt_patricia EXPORT dt_patricia-targets)
install(EXPORT dt_patricia-targets
    FILE dt_patricia-targets.cmake
    NAMESPACE dt_patricia::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dt_patricia
)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dt_patricia-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/dt_patricia-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dt_patricia
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/dt_patricia-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/dt_patricia-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/dt_patricia-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dt_patricia
)
