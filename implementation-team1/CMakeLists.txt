cmake_minimum_required(VERSION 2.8.4)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Release, Debug" FORCE)
endif ()


project(mlalign)

# Compiler Warnings

if (CMAKE_CXX_COMPILER_ID STREQUAL GNU OR
    CMAKE_CXX_COMPILER_ID STREQUAL Clang)

    set(specific_flags "-std=c++11 -Wall  -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimport  -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized") 

elseif (CMAKE_CXX_COMPILER_ID STREQUAL Clang)

    set(specific_flags "-std=c++11 -Weverything")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL MSVC)

    set(specific_flags "/W4 /WX /EHsc")

endif ()


# Include dirs

# This is so that libraries can include this as a sub project
set(MLALIGN_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
get_directory_property(hasParent PARENT_DIRECTORY)
if (hasParent)

    # Exports the variable to the parent scope
    set(MLALIGN_INCLUDE_DIRS ${MLALIGN_INCLUDE_DIRS} PARENT_SCOPE)

endif ()


if (NOT SIMD_ARCH)
  # Try some auto-detection magic
  include(SimdMultiarch.cmake)
  simdpp_get_runnable_archs(NATIVE_ARCHS)
  list(GET NATIVE_ARCHS -1 SIMD_ARCH) # The last item is the highest vectorizing architecture
endif ()

if (SIMD_ARCH MATCHES "X86_AVX")
  add_definitions("-mavx")
elseif (SIMD_ARCH MATCHES "X86_SSE3")
  add_definitions("-msse3")
elseif (SIMD_ARCH MATCHES "X86_Scalar")
  # This is just the fallback behavior
else ()
  message(ERROR "Unknown architecture: SIMD_ARCH=" ${SIMD_ARCH})
endif ()

message(STATUS "Compiling for " ${SIMD_ARCH})
add_definitions("-DSIMD_ARCH_${SIMD_ARCH}")


# Lib configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${specific_flags}")
set(LIB_SOURCE_FILES src/mlalign.cpp)
include_directories(${MLALIGN_INCLUDE_DIRS} ${Vc_INCLUDE_DIR})
add_library(libmlalign ${LIB_SOURCE_FILES})
target_link_libraries(libmlalign ${Vc_LIBRARIES})


# Console wrapper configuration

set(EXE_SOURCE_FILES src/main.cpp)
add_executable(mlalign ${EXE_SOURCE_FILES})
target_link_libraries(mlalign libmlalign)
