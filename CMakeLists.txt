cmake_minimum_required(VERSION 2.8.12)

project(mindex)

set(MAX_KMER_SIZE "32")
set(MAX_GMER_SIZE "32")
SET(CMAKE_C_FLAGS "-std=c99 -DXXH_NAMESPACE=MINDEX_HASH_")
SET(CMAKE_CXX_FLAGS "-std=c++11 -DXXH_NAMESPACE=MINDEX_HASH_")

set_property(SOURCE xxhash.c APPEND_STRING PROPERTY COMPILE_FLAGS " -Wall -Wextra -Wcast-qual -Wcast-align -Wshadow -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement -Wstrict-prototypes -Wundef")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("debug mode")
	add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_BUILD_TYPE MATCHES Profile)
		message("profile mode")
		add_compile_options(-pg)
		set(CMAKE_SHARED_LINKER_FLAGS "-pg")
		set(CMAKE_EXE_LINKER_FLAGS "-pg")
	else(CMAKE_BUILD_TYPE MATCHES Profile)
		message("release mode")
		add_compile_options(-O3)
	endif(CMAKE_BUILD_TYPE MATCHES Profile)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
