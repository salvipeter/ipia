# -*- mode: cmake -*-

cmake_minimum_required(VERSION 3.1.0)

project(ipia)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(LIBGEOM_DIR "../libgeom") # https://github.com/salvipeter/libgeom
include_directories("${LIBGEOM_DIR}")
link_directories("${LIBGEOM_DIR}/release")

add_library(ipia STATIC ipia.cc)

target_link_libraries(ipia "geom")

set(LIBDC_DIR "../dual-contouring") # https://github.com/salvipeter/dual-contouring
include_directories("${LIBGEOM_DIR}" "${LIBDC_DIR}")
link_directories("${LIBGEOM_DIR}/release" "${LIBDC_DIR}/build")

add_executable(ipia-test ipia-test.cc)

target_link_libraries(ipia-test ipia "dualcontour")
