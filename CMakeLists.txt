cmake_minimum_required(VERSION 2.8)

project(libsaxs C)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

# Version
set(saxs_version_major 1)
set(saxs_version_minor 0)
set(saxs_version_patch 0)
set(saxs_version
  "${saxs_version_major}.${saxs_version_minor}.${saxs_version_patch}")

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always HEAD
    COMMAND sed -e "s/-/./g"
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GIT_VERSION "git")
endif()

# build options
option( PROFILE	                        "Generate profiling information"                OFF )
option( OPEN_MP                         "Build With OpenMP"                             ON )

set( CMAKE_INSTALL_PREFIX "$ENV{HOME}"                      CACHE PATH "Base directory for installation" FORCE )
set( LIB_INSTALL_DIR      "${CMAKE_INSTALL_PREFIX}/lib"     CACHE PATH "Where to install libraries to." )
set( INCLUDE_INSTALL_DIR  "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Where to install headers to." )

if(PROFILE)
  set(CMAKE_C_FLAGS "-pg ${CMAKE_C_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -fPIC")
if(OPEN_MP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
endif()

include(CheckCCompilerFlag)
check_c_compiler_flag(-flto HAS_LINK_TIME_OPTIMIZATION)
if(HAS_LINK_TIME_OPTIMIZATION)
  set(CMAKE_C_FLAGS_RELEASE
    "${CMAKE_C_FLAGS_RELEASE} -flto")
endif()

# Libraries
find_package(libmol REQUIRED)
include_directories(${LIBMOL_INCLUDE_DIRS})

# Targets
add_library(saxs
  src/saxs_profile.c
  src/saxs_atom.c
  src/form_factor_table.c
  src/radial_distribution.c
  src/rgyration.c
  src/saxs_score.c
  src/saxs_utils.c)
target_link_libraries(saxs
  ${LIBMOL_LIBRARIES}
  m z)

install(TARGETS saxs
  LIBRARY DESTINATION ${LIB_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
file(GLOB header_files "src/*.h")
install(FILES ${header_files} DESTINATION ${INCLUDE_INSTALL_DIR}/saxs.${saxs_version})

# Subdirectories
add_subdirectory(tools)
