#include "geometry_tests.h"
#include "test_assert.h"

#include "../jtk/geometry.h"

namespace jtk
  {
  
  void test_stitch()
    {
    std::vector<jtk::vec3<float>> vertices;
    vertices.emplace_back(0, 0, 0);
    vertices.emplace_back(1, 0, 0);
    vertices.emplace_back(1, 1, 0);
    vertices.emplace_back(0, 1, 0);
    vertices.emplace_back(0.5, 0.5, 0);
    vertices.emplace_back(0.51, 0.5, 0);
    
    std::vector<jtk::vec3<uint32_t>> triangles;
    triangles.emplace_back(0, 1, 4);
    triangles.emplace_back(1, 2, 5);
    triangles.emplace_back(2, 5, 3);
    triangles.emplace_back(3, 4, 0);
    
    TEST_ASSERT(!stitch_points(vertices, triangles, 0.0001f));
    TEST_ASSERT(stitch_points(vertices, triangles, 0.1f));
    TEST_EQ(5, vertices.size());
    }
  
  }

void run_all_geometry_tests()
  {

  }
