
file(GLOB RANDOMWALK_SOURCES "${SOURCES_DIR}/randomwalk/*.cpp")
file(GLOB RANDOMWALK_HEADERS "${SOURCES_DIR}/randomwalk/*.hpp")

set(RANDOMWALK_FILES ${RANDOMWALK_SOURCES} ${RANDOMWALK_HEADERS})

set(ALL_FILES
  ${RANDOMWALK_FILES}
)

find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(randomwalk ${ALL_FILES})
target_link_libraries(randomwalk PUBLIC Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(randomwalk PUBLIC . Eigen3::Eigen)