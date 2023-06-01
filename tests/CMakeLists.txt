find_package(Catch2 REQUIRED)

add_executable(test_pw ${TESTS_DIR}/test.cpp)

target_link_libraries(test_pw PRIVATE randomwalk Catch2::Catch2WithMain)