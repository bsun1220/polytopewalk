find_package(Catch2 REQUIRED)

add_executable(test_pw ${TESTS_DIR}/test.cpp)
add_executable(test_sparse ${TESTS_DIR}/testsparse.cpp)

target_link_libraries(test_pw PRIVATE randomwalk Catch2::Catch2WithMain)
target_link_libraries(test_sparse PRIVATE sparserandomwalk Catch2::Catch2WithMain)