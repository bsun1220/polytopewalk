
file(GLOB DENSE_SOURCES "${SOURCES_DIR}/dense/*.cpp")
file(GLOB DENSE_HEADERS "${SOURCES_DIR}/dense/*.hpp")

file(GLOB SPARSE_SOURCES "${SOURCES_DIR}/sparse/*.cpp")
file(GLOB SPARSE_HEADERS "${SOURCES_DIR}/sparse/*.hpp")

file(GLOB UTILS_SOURCES "${SOURCES_DIR}/utils/*.cpp")
file(GLOB UTILS_HEADERS "${SOURCES_DIR}/utils/*.hpp")

set(DENSE_FILES ${DENSE_SOURCES} ${DENSE_HEADERS})
set(SPARSE_FILES ${SPARSE_SOURCES} ${SPARSE_HEADERS})
set(UTILS_FILES ${UTILS_SOURCES} ${UTILS_HEADERS})

set(ALL_FILES
  ${DENSE_FILES}
  ${SPARSE_FILES}
  ${UTILS_FILES}
)

find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(dense ${DENSE_FILES})
target_link_libraries(dense PUBLIC Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(dense PUBLIC . Eigen3::Eigen)

add_library(sparse ${SPARSE_FILES})
target_link_libraries(sparse PUBLIC Eigen3::Eigen)
target_include_directories(sparse PUBLIC . Eigen3::Eigen)

add_library(utils ${UTILS_FILES})
target_link_libraries(utils PUBLIC Eigen3::Eigen ifopt::ifopt_ipopt dense sparse)
target_include_directories(utils PUBLIC . Eigen3::Eigen)