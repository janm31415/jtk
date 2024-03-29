

if (${JTK_THREADING} STREQUAL "tbb")
add_definitions(-D_ENABLE_TBB)
endif (${JTK_THREADING} STREQUAL "tbb")

if (${JTK_THREADING} STREQUAL "ppl")
add_definitions(-D_ENABLE_PPL)
endif (${JTK_THREADING} STREQUAL "ppl")

if (${JTK_THREADING} STREQUAL "std")
add_definitions(-D_ENABLE_THREADS)
endif (${JTK_THREADING} STREQUAL "std")

if (JTK_MAT_PARALLEL)
add_definitions(-D_JTK_MAT_PARALLEL)
endif (JTK_MAT_PARALLEL)

if (JTK_NO_SIMD)
add_definitions(-D_JTK_NO_SIMD)
endif (JTK_NO_SIMD)

if (JTK_MAT_SIMD)
add_definitions(-D_JTK_MAT_SIMD)
endif (JTK_MAT_SIMD)


if (${JTK_TARGET} STREQUAL "arm")
add_definitions(-D_JTK_FOR_ARM)
endif (${JTK_TARGET} STREQUAL "arm")
