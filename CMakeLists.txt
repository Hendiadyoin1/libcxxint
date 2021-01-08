cmake_minimum_required (VERSION 2.6)
project (cint C CXX)

message("TOOLCHAINFILE: ${CMAKE_TOOLCHAIN_FILE}")

set(cint_VERSION_MAJOR "4")
set(cint_VERSION_MINOR "0")
set(cint_VERSION_PATCH "8")
set(cint_VERSION_TWEAK "0")
set(cint_VERSION "${cint_VERSION_MAJOR}.${cint_VERSION_MINOR}.${cint_VERSION_PATCH}")
set(cint_SOVERSION "${cint_VERSION_MAJOR}")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -DNDEBUG")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()
set(CMAKE_VERBOSE_MAKEFILE OFF)

if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} VERSION_GREATER 2.8.3)
  if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION} VERSION_GREATER 3.0.0)
    cmake_policy(SET CMP0042 NEW) 
  endif()
  include(GNUInstallDirs)
else ()
  set(CMAKE_INSTALL_LIBDIR lib CACHE PATH "lib")
  set(CMAKE_INSTALL_INCLUDEDIR include CACHE PATH "include")
endif ()
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(QUADMATH)
if(QUADMATH_FOUND)
  set(HAVE_QUADMATH_H 1 CACHE INTERNAL "HAVE QUADMATH")
endif()



include(CheckFunctionExists)
set(CMAKE_REQUIRED_INCLUDES math.h )
set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(expl HAVE_EXPL)
check_function_exists(sqrtl HAVE_SQRTL)
#set(CMAKE_REQUIRED_INCLUDES quadmath.h)
#set(CMAKE_REQUIRED_LIBRARIES quadmath)
#check_function_exists(fabsq HAVE_QUADMATH_H)

include_directories(
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/config.hpp.in"
  "${PROJECT_BINARY_DIR}/src/config.hpp")

configure_file(
  "${PROJECT_SOURCE_DIR}/include/cint.hpp.in"
  "${PROJECT_BINARY_DIR}/include/cint.hpp")

set(cintSrc 
  src/c2f.cpp src/cart2sph.cpp src/cint1e.cpp src/cint2e.cpp src/cint_bas.cpp
  src/fblas.cpp src/g1e.cpp src/g2e.cpp src/misc.cpp src/optimizer.cpp
  src/rys_roots.cpp
  src/cint2c2e.cpp src/g2c2e.cpp src/cint3c2e.cpp src/g3c2e.cpp
  src/cint3c1e.cpp src/g3c1e.cpp src/breit.cpp src/cint4c1e.cpp src/g4c1e.cpp
  src/cint1e_a.cpp src/cint3c1e_a.cpp
  src/autocode/breit1.cpp src/autocode/dkb.cpp src/autocode/gaunt1.cpp
  src/autocode/grad1.cpp src/autocode/grad2.cpp src/autocode/hess.cpp
  src/autocode/int3c1e.cpp src/autocode/int3c2e.cpp src/autocode/intor1.cpp
  src/autocode/intor2.cpp src/autocode/intor3.cpp src/autocode/intor4.cpp
  src/autocode/deriv3.cpp)

#option(WITH_RANGE_COULOMB "Range separated Coulomb" on)
if(WITH_RANGE_COULOMB)
# defined in config.h
#  add_definitions(-DWITH_RANGE_COULOMB)
  message("Enabled WITH_RANGE_COULOMB")
  set(cintSrc ${cintSrc} src/erfc_roots.cpp src/polyfits.cpp)
endif(WITH_RANGE_COULOMB)

if(WITH_COULOMB_ERF)
  set(cintSrc ${cintSrc} src/g2e_coulerf.cpp src/cint2e_coulerf.cpp)
  add_definitions(-DWITH_COULOMB_ERF)
  message("Enabled WITH_COULOMB_ERF")
endif(WITH_COULOMB_ERF)

if(WITH_F12)
  set(cintSrc ${cintSrc} src/g2e_f12.cpp src/stg_roots.cpp src/cint2e_f12.cpp)
  add_definitions(-DWITH_F12)
  message("Enabled WITH_F12")
endif(WITH_F12)

if(WITH_GTG)
  set(cintSrc ${cintSrc} src/g2e_gtg.cpp src/cint2e_gtg.cpp src/cint3c2e_gtg.cpp
    src/cint2c2e_gtg.cpp)
  add_definitions(-DWITH_GTG)
  message("Enabled WITH_GTG")
endif(WITH_GTG)

if(PYPZPX)
  add_definitions(-DPYPZPX)
  message("P orbitals convention (py, pz, px)")
else()
  message("P orbitals convention (px, py, pz)")
endif(PYPZPX)

if(MIN_EXPCUTOFF)
  message("Minimal cutoff in exponent for integral screening ${MIN_EXPCUTOFF}")
  add_definitions(-DMIN_EXPCUTOFF=${MIN_EXPCUTOFF})
endif()

if(KEEP_GOING)
  message("Whether to trigger hard exit for numerical issues in Rys quadrature")
  add_definitions(-DKEEP_GOING)
endif()

option(WITH_FORTRAN "Fortran interface" on)
if(WITH_FORTRAN)
  add_definitions(-DWITH_FORTRAN)
  message("Enable Fortran interface")
else()
  message("Exclude Fortran interface")
endif(WITH_FORTRAN)

option(WITH_CINT2_INTERFACE "Old libcint (version 2) interface" on)
if(WITH_CINT2_INTERFACE)
  add_definitions(-DWITH_CINT2_INTERFACE)
  message("Enable old cint (version 2) interface")
else()
  message("Exclude old cint (version 2) interface")
endif(WITH_CINT2_INTERFACE)

option(BUILD_SHARED_LIBS "build shared libraries" 1)
option(ENABLE_EXAMPLE "build examples" 0)
option(ENABLE_TEST "build tests" 1)
option(ENABLE_STATIC "Enforce static library build" 1)
if(QUICK_TEST)
  set(RUN_QUICK_TEST --quick)
endif()


if(ENABLE_STATIC)
  set(BUILD_SHARED_LIBS 0)
endif()

# if(BUILD_SHARED_LIBS)
#   find_package(OpenBLAS CONFIG REQUIRED COMPONENTS shared)
# else()
#   find_package(OpenBLAS CONFIG REQUIRED COMPONENTS static)
# endif()
find_package(BLAS REQUIRED)

add_library(cint ${cintSrc})
set_target_properties(cint PROPERTIES
  VERSION ${cint_VERSION}
  SOVERSION ${cint_SOVERSION}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

target_link_libraries(cint "-lm")
message("using blas: ${BLAS_LIBRARIES}")
target_link_libraries(cint ${BLAS_LIBRARIES})

if(QUADMATH_FOUND)
  target_link_libraries(cint quadmath)
endif()

set(CintHeaders
  ${PROJECT_SOURCE_DIR}/include/cint_funcs.hpp
  ${PROJECT_BINARY_DIR}/include/cint.hpp)

install(TARGETS cint DESTINATION "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" COMPONENT "lib")
install(FILES ${CintHeaders} DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} COMPONENT "dev")


if(ENABLE_EXAMPLE)
  enable_language(Fortran)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(HAVE_OPENMP 1)
  else ()
    set(OpenMP_CXX_FLAGS " ")
  endif()
  add_subdirectory(examples)
endif()

if(BUILD_SHARED_LIBS AND ENABLE_TEST)
  find_package(PythonInterp)
  message(STATUS "Found python  ${PYTHON_EXECUTABLE}")
  if(PYTHONINTERP_FOUND)
    enable_testing()
    add_test(NAME cinttest
      COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/testsuite/test_cint.py ${RUN_QUICK_TEST})
    add_test(NAME cint3c2etest
      COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/testsuite/test_3c2e.py ${RUN_QUICK_TEST})
  endif()
endif()

include(CPackComponent)

set(CPACK_PACKAGE_CHECKSUM "SHA512")
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_COMPRESSION_TYPE "xz")
set(CPACK_PACKAGE_VENDOR "Qiming @sunqm Sun")
set(CPACK_PACKAGE_CONTACT "https://github.com/sunqm/libcint")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/sunqm/libcint")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/doc/libcint.CPack.txt")

set(CPACK_PACKAGE_VERSION_MAJOR "${cint_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${cint_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${cint_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_TWEAK "${cint_VERSION_TWEAK}")

set(CPACK_COMPONENTS_GROUPING "IGNORE")

set(CPACK_MONOLITHIC_INSTALL OFF)
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_RPM_COMPONENT_INSTALL ON)
set(CPACK_NSIS_COMPONENT_INSTALL ON)

set("CPACK_DEBIAN_LIB_PACKAGE_NAME" "libcint")
set("CPACK_DEBIAN_DEV_PACKAGE_NAME" "${CPACK_DEBIAN_LIB_PACKAGE_NAME}-dev")

set(CPACK_DEBIAN_LIB_PACKAGE_DEPENDS "libopenblas0-pthread")
set(CPACK_DEBIAN_DEV_PACKAGE_DEPENDS "${CPACK_DEBIAN_LIB_PACKAGE_NAME},libopenblas-pthread-dev")
if(HAVE_QUADMATH)
  set(CPACK_DEBIAN_LIB_PACKAGE_DEPENDS "${CPACK_DEBIAN_LIB_PACKAGE_DEPENDS},libquadmath0")
endif()

set(CPACK_COMPONENTS_GROUPING "IGNORE")

include(CPack)
