set(JTK_THREADING "std" CACHE STRING  "Choose your threading library: tbb is Intel's Threading Building Blocks library, ppl is Windows' concurrency library, std uses std::thread, and none does not use threading.")
 
set_property(CACHE JTK_THREADING PROPERTY STRINGS tbb ppl std none)

option(JTK_MAT_PARALLEL "Certain matrix operations (e.g. conjugate gradients) run in parallel" ON)

option(JTK_NO_SIMD "No SIMD instructions" OFF)

option(JTK_MAT_SIMD "Run matrix operations with SIMD" OFF)

set(JTK_TARGET "x64" CACHE STRING  "Choose your target platform: x64 compatible pc or arm.")
set_property(CACHE JTK_TARGET PROPERTY STRINGS x64 arm)

set(HDRS
alsc.h
clamp.h
concurrency.h
containers.h
deformation.h
file_utils.h
fitting.h
geometry.h
halffloat.h
icp.h
image.h
log.h
octree.h
opengl.h
mat.h
pipe.h
ply.h
point_tree.h
qbvh.h
rand.h
render.h
sse2neon.h
timer.h
utf8.h
vec.h
window.h
)
	
set(SRCS
jtk.cpp
)

if (UNIX)
set(CMAKE_C_FLAGS_DEBUG "-msse4.1")
set(CMAKE_CXX_FLAGS_DEBUG "-msse4.1")
set(CMAKE_C_FLAGS_RELEASE "-msse4.1")
set(CMAKE_CXX_FLAGS_RELEASE "-msse4.1")
endif (UNIX)

if (APPLE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
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
