cmake_minimum_required(VERSION 2.8)
project(Couette)

IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release)
ENDIF()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(
${PROJECT_SOURCE_DIR}/include
)

aux_source_directory(src/ DIR_SRCS)

add_executable(Couette ${DIR_SRCS})
