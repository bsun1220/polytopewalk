#find_package(Catch2 REQUIRED)
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

add_executable(test ${TESTS_DIR}/test.cpp)
target_link_libraries(test PRIVATE utils dense sparse Catch2::Catch2WithMain)