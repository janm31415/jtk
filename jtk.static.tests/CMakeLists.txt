if (WIN32)
  set(TBB_INCLUDE_DIR "C:/Program Files/TBB/include")
  set(TBB_LIBRARIES "C:/Program Files/TBB/lib/intel64/vc14/tbb.lib" "C:/Program Files/TBB/lib/intel64/vc14/tbb_debug.lib")
endif (WIN32)

if (UNIX)
  set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include")
  set(TBB_LIBRARIES "/usr/lib/x86_64-linux-gnu/libtbb.so" "/usr/lib/x86_64-linux-gnu/libtbb.so.2")  
endif (UNIX)

if (APPLE)
  set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include")
  set(TBB_LIBRARIES "/usr/local/lib/libtbb.dylib")
endif (APPLE)

set(HDRS
    )
	
set(SRCS
static_tests.cpp
static_tests_2.cpp
)

if (WIN32)
set(CMAKE_C_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_CXX_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_C_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /bigobj")

set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS_DEBUG} /WX")
set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS_RELEASE} /WX")
endif (WIN32)

if (UNIX)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pthread -mavx2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pthread -mavx2")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pthread -mavx2")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pthread -mavx2")
endif (UNIX)

if (APPLE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
endif (APPLE)

# general build definitions
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
if (WIN32)
add_definitions(-D_UNICODE)
add_definitions(-DUNICODE)
endif (WIN32)
include ("../jtk/jtk.cmake")

add_definitions(-DMEMORY_LEAK_TRACKING)

add_library(jtk.static.tests STATIC ${SRCS})
source_group("Source Files" FILES ${srcs})	

if (UNIX)

set(XLIBINCLUDE "")
set(XLIBLIBRARY "X11.so")

if (APPLE)
set(XLIBINCLUDE "/usr/X11/include/")
set(XLIBLIBRARY "/usr/X11/lib/libX11.dylib")
endif (APPLE)

target_include_directories(jtk.static.tests
   PRIVATE
   ${XLIBINCLUDE}
   )
   
target_link_libraries(jtk.static.tests
    PRIVATE
    ${XLIBLIBRARY}
    )
endif (UNIX)

if (${JTK_THREADING} STREQUAL "tbb")
  target_include_directories(jtk.static.tests
      PRIVATE
      ${TBB_INCLUDE_DIR}
      )
  target_link_libraries(jtk.static.tests
      PRIVATE
      ${TBB_LIBRARIES}
      )
endif (${JTK_THREADING} STREQUAL "tbb")
