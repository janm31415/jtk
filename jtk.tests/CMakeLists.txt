if (WIN32)
  set(TBB_INCLUDE_DIR "C:/Program Files/TBB/include")
  set(TBB_LIBRARIES "C:/Program Files/TBB/lib/intel64/vc14/tbb.lib" "C:/Program Files/TBB/lib/intel64/vc14/tbb_debug.lib")
endif (WIN32)

if (UNIX)
  set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include")
  set(TBB_LIBRARIES "/usr/lib/x86_64-linux-gnu/libtbb.so" "/usr/lib/x86_64-linux-gnu/libtbb.so.2")  
endif (UNIX)

set(HDRS
file_utils_tests.h
mat_tests.h
qbvh_tests.h
test_assert.h
vec_tests.h
    )
	
set(SRCS
file_utils_tests.cpp
mat_tests.cpp
qbvh_tests.cpp
test_assert.cpp
test.cpp
vec_tests.cpp
)

if (WIN32)
set(CMAKE_C_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi /bigobj")
set(CMAKE_CXX_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi /bigobj")
set(CMAKE_C_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG /bigobj")
set(CMAKE_CXX_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG /bigobj")
endif (WIN32)

if (UNIX)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -msse4.1")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -msse4.1")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -msse4.1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse4.1")
endif (UNIX)

# general build definitions
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_definitions(-D_MAT_PARALLEL)
#add_definitions(-D_ENABLE_TBB)
#add_definitions(-D_ENABLE_PPL)
add_definitions(-D_ENABLE_THREADS)

add_definitions(-DMEMORY_LEAK_TRACKING)

add_executable(jtk.tests ${HDRS} ${SRCS})
source_group("Header Files" FILES ${hdrs})
source_group("Source Files" FILES ${srcs})

#if (WIN32)
#add_custom_command(TARGET jtk.tests POST_BUILD 
#   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/data" ${CMAKE_CURRENT_BINARY_DIR}/data)
#endif (WIN32)

#if (UNIX)
#add_custom_command(TARGET jtk.tests POST_BUILD 
#   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/data" "$<TARGET_FILE_DIR:libskiwi.tests>/data")
#endif (UNIX)



target_include_directories(jtk.tests
    PRIVATE
    ${TBB_INCLUDE_DIR}
    )
	
target_link_libraries(jtk.tests
    PRIVATE
    ${TBB_LIBRARIES}
    )	