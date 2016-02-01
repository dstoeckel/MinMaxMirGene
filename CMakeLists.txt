project(MinMaxMirGene)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(SOURCES
	main.cpp
	MinMaxProblem.cpp
	TargetMappings.cpp
)

set(HEADERS
	MinMaxProblem.h
	TargetMappings.h
)

include_directories(
	${CPLEX_INCLUDE_DIR}
)

link_directories(
	${CPLEX_LIBRARY_DIR}
)

set(COMPILER_FLAGS
	"-Wall -std=c++14 -pedantic -fpie"
)

set(LINK_FLAGS
	"-Wl,--as-needed -pie"
)

add_executable(minMaxMirGene ${SOURCES} ${HEADERS})
target_link_libraries(minMaxMirGene cplex1262 dl)
set_target_properties(minMaxMirGene PROPERTIES
	COMPILE_FLAGS ${COMPILER_FLAGS}
	LINK_FLAGS ${LINK_FLAGS}
)

install(
	TARGETS minMaxMirGene
	RUNTIME DESTINATION bin
)

