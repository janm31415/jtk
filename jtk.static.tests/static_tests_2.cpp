#pragma warning(push)
#pragma warning(disable:4505)

#define JTK_IMAGE_IMPLEMENTATION
#define JTK_IMAGE_STATIC
#include "../jtk/image.h"

//#define JTK_DEFORMATION_IMPLEMENTATION
//#define JTK_DEFORMATION_STATIC
//#include "../jtk/deformation.h"

#define JTK_GEOMETRY_IMPLEMENTATION
#define JTK_GEOMETRY_STATIC
#include "../jtk/geometry.h"

#define JTK_ALSC_IMPLEMENTATION
#define JTK_ALSC_STATIC
#include "../jtk/alsc.h"

#define JTK_FILE_UTILS_IMPLEMENTATION
#define JTK_FILE_UTILS_STATIC
#include "../jtk/file_utils.h"

#define JTK_ICP_IMPLEMENTATION
#define JTK_ICP_STATIC
#include "../jtk/icp.h"

#define JTK_PIPE_IMPLEMENTATION
#define JTK_PIPE_STATIC
#include "../jtk/pipe.h"

#define JTK_PLY_IMPLEMENTATION
#define JTK_PLY_STATIC
#include "../jtk/ply.h"

#define JTK_QBVH_IMPLEMENTATION
#define JTK_QBVH_STATIC
#include "../jtk/qbvh.h"

#ifdef _WIN32
#define JTK_WINDOW_IMPLEMENTATION
#define JTK_WINDOW_STATIC
#include "../jtk/window.h"
#endif

#define JTK_HALFFLOAT_IMPLEMENTATION
#define JTK_HALFFLOAT_STATIC
#include "../jtk/halffloat.h"

#pragma warning(pop)
