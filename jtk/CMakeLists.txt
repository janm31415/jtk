
set(HDRS
alsc.h
clamp.h
concurrency.h
file_utils.h
fitting.h
geometry.h
image.h
mat.h
pipe.h
qbvh.h
render.h
sse2neon.h
timer.h
utf8.h
vec.h
)
	
set(SRCS
jtk.cpp
)

if (WIN32)
set(CMAKE_C_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_CXX_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_C_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
endif (WIN32)

if (UNIX)
set(CMAKE_C_FLAGS_DEBUG "-msse4.1")
set(CMAKE_CXX_FLAGS_DEBUG "-msse4.1")
set(CMAKE_C_FLAGS_RELEASE "-msse4.1")
set(CMAKE_CXX_FLAGS_RELEASE "-msse4.1")
endif (UNIX)

if (APPLE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
endif (APPLE)

# general build definitions
add_definitions(-D_SCL_SECURE_NO_WARNINGS)

add_library(jtk STATIC ${HDRS} ${SRCS})
source_group("Header Files" FILES ${HDRS})
source_group("Source Files" FILES ${SRCS})

target_include_directories(jtk
    PRIVATE
    )	
	
target_link_libraries(jtk
    PRIVATE	
    )	
