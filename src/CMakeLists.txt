
file(GLOB RANDOMWALK_SOURCES "${SOURCES_DIR}/randomwalk/*.cpp")
file(GLOB RANDOMWALK_HEADERS "${SOURCES_DIR}/randomwalk/*.hpp")

file(GLOB SPARSERANDOMWALK_SOURCES "${SOURCES_DIR}/sparserandomwalk/*.cpp")
file(GLOB SPARSERANDOMWALK_HEADERS "${SOURCES_DIR}/sparserandomwalk/*.hpp")

file(GLOB UTILS_SOURCES "${SOURCES_DIR}/utils/*.cpp")
file(GLOB UTILS_HEADERS "${SOURCES_DIR}/utils/*.hpp")

file(GLOB CONSTRAINTWALK_SOURCES "${SOURCES_DIR}/constraintwalk/*.cpp")
file(GLOB CONSTRAINTWALK_HEADERS "${SOURCES_DIR}/constraintwalk/*.hpp")

set(RANDOMWALK_FILES ${RANDOMWALK_SOURCES} ${RANDOMWALK_HEADERS})
set(SPARSERANDOMWALK_FILES ${SPARSERANDOMWALK_SOURCES} ${SPARSERANDOMWALK_HEADERS})
set(UTILS_FILES ${UTILS_SOURCES} ${UTILS_HEADERS})
set(CONSTRAINTWALK_FILES ${CONSTRAINTWALK_SOURCES} ${CONSTRAINTWALK_HEADERS})

set(ALL_FILES
  ${RANDOMWALK_FILES}
  ${SPARSERANDOMWALK_FILES}
  ${UTILS_FILES}
  ${CONSTRAINTWALK_FILES}
)

find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)

#add_library(randomwalk ${RANDOMWALK_FILES})
#target_link_libraries(randomwalk PUBLIC Eigen3::Eigen ifopt::ifopt_ipopt)
#target_include_directories(randomwalk PUBLIC . Eigen3::Eigen)

#add_library(sparserandomwalk ${SPARSERANDOMWALK_FILES})
#target_link_libraries(sparserandomwalk PUBLIC Eigen3::Eigen)
#target_include_directories(sparserandomwalk PUBLIC . Eigen3::Eigen)

add_library(utils ${UTILS_FILES})
target_link_libraries(utils PUBLIC Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(utils PUBLIC . Eigen3::Eigen)

add_library(constraintwalk ${CONSTRAINTWALK_FILES})
target_link_libraries(constraintwalk PUBLIC Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(constraintwalk PUBLIC . Eigen3::Eigen)