file(GLOB sources *.cpp *c)
file(GLOB headers *.hpp *.h)

list(REMOVE_ITEM sources mindex_main.cpp)

add_definitions( -DMAX_KMER_SIZE=${MAX_KMER_SIZE} -DMAX_GMER_SIZE=${MAX_GMER_SIZE})

add_library(mindex_core ${sources} ${headers})
target_include_directories(mindex_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mindex mindex_main.cpp)

find_package( Threads REQUIRED )
target_link_libraries(mindex mindex_core pthread)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)

    if (UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(mindex PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)

find_package( ZLIB REQUIRED )
target_link_libraries( mindex ${ZLIB_LIBRARIES})

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
        target_link_libraries(mindex librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        target_link_libraries(mindex rt)
    endif()
endif(LINK MATCHES static)

install(TARGETS mindex DESTINATION bin)
