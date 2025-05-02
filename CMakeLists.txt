cmake_minimum_required(VERSION 3.10)

project(Algebra LANGUAGES C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR}/inc)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/main.c
    ${PROJECT_SOURCE_DIR}/src/algebra.c
)

add_executable(algebra ${SOURCES})