cmake_minimum_required(VERSION 3.20)

project(rtfed3d LANGUAGES C CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(RTFED3D_BUILD_EXAMPLES "Build example simulation drivers" ON)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(PkgConfig QUIET)

set(RTFED3D_HDF5_LIBRARIES)
set(RTFED3D_HDF5_INCLUDE_DIRS)

if(PkgConfig_FOUND)
  pkg_search_module(HDF5 REQUIRED IMPORTED_TARGET hdf5-openmpi hdf5-mpi)
  set(RTFED3D_HDF5_LIBRARIES PkgConfig::HDF5)
  set(RTFED3D_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
else()
  find_package(HDF5 REQUIRED COMPONENTS C)
  if(NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Parallel HDF5 is required, but the discovered HDF5 build is not MPI-enabled.")
  endif()
  set(RTFED3D_HDF5_LIBRARIES ${HDF5_LIBRARIES})
  set(RTFED3D_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
endif()

find_path(BLITZ_INCLUDE_DIR NAMES blitz/array.h)
find_library(BLITZ_LIBRARY NAMES blitz)

if(NOT BLITZ_INCLUDE_DIR OR NOT BLITZ_LIBRARY)
  message(FATAL_ERROR "Blitz++ was not found. Set BLITZ_INCLUDE_DIR and BLITZ_LIBRARY, or install Blitz++ in a standard location.")
endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/base/common.hpp")
  message(FATAL_ERROR "Vendored support code is missing from base/. Ensure the repository contains the base/ directory.")
endif()

add_library(blitz::blitz UNKNOWN IMPORTED)
set_target_properties(blitz::blitz PROPERTIES
  IMPORTED_LOCATION "${BLITZ_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${BLITZ_INCLUDE_DIR}"
)

set(RTFED3D_CORE_SOURCES
  src/bc_periodic.cpp
  src/global.cpp
  src/hdf5util.cpp
  src/integrator.cpp
  src/io_hdf5.cpp
  src/mc2.cpp
  src/solver.cpp
)

add_library(rtfed3d_core STATIC ${RTFED3D_CORE_SOURCES})
target_include_directories(rtfed3d_core
  PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/base
    ${RTFED3D_HDF5_INCLUDE_DIRS}
)
target_link_libraries(rtfed3d_core
  PUBLIC
    MPI::MPI_CXX
    blitz::blitz
    ${RTFED3D_HDF5_LIBRARIES}
)

if(RTFED3D_BUILD_EXAMPLES)
  set(RTFED3D_EXAMPLE_NAMES
    blastwave
    circular
    diffusion
    reconnection
    riemann1d
    sheets
    tearing
    vortex
  )

  foreach(example_name IN LISTS RTFED3D_EXAMPLE_NAMES)
    add_executable(${example_name} examples/${example_name}.cpp)
    target_link_libraries(${example_name} PRIVATE rtfed3d_core)
    set_target_properties(${example_name} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples"
    )
    add_custom_command(TARGET ${example_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/examples/${example_name}.cfg"
        "$<TARGET_FILE_DIR:${example_name}>/${example_name}.cfg"
      VERBATIM
    )
  endforeach()
endif()

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)

  add_test(NAME riemann1d_help
    COMMAND $<TARGET_FILE:riemann1d> --help
  )
  set_tests_properties(riemann1d_help PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples
  )

  add_test(NAME riemann1d_small_regression
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/run_and_compare_hdf5.py
      --executable $<TARGET_FILE:riemann1d>
      --config ${CMAKE_SOURCE_DIR}/examples/riemann1d_small.cfg
      --reference ${CMAKE_SOURCE_DIR}/tests/data/riemann1d_small_reference.h5
      --output-name riemann1d_small.h5
  )
endif()
