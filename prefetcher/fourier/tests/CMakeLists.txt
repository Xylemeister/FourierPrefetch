cmake_minimum_required(VERSION 3.14)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

FetchContent_MakeAvailable(googletest)

enable_testing()

#   add_library(lib
#  ../transform.cpp
#  )

# --- Test executable ---
add_executable(tests
    test_transform.cpp
)

target_include_directories(tests PRIVATE ..)


target_link_libraries(tests
   # lib
    GTest::gtest_main
)

# Auto-discover and register tests with CTest
include(GoogleTest)
gtest_discover_tests(tests)
