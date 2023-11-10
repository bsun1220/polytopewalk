#find_package(Catch2 REQUIRED)
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

add_executable(test_pw ${TESTS_DIR}/test.cpp)
add_executable(test_sparse ${TESTS_DIR}/testsparse.cpp)

target_link_libraries(test_pw PRIVATE randomwalk utils Catch2::Catch2WithMain)
target_link_libraries(test_sparse PRIVATE sparserandomwalk utils Catch2::Catch2WithMain)