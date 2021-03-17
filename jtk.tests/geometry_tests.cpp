#include "geometry_tests.h"
#include "test_assert.h"

#define JTK_GEOMETRY_IMPLEMENTATION
#include "../jtk/geometry.h"

namespace jtk
  {
  
  void test_stitch()
    {
    std::vector<jtk::vec3<float>> vertices;
    vertices.emplace_back(0.f, 0.f, 0.f);
    vertices.emplace_back(1.f, 0.f, 0.f);
    vertices.emplace_back(1.f, 1.f, 0.f);
    vertices.emplace_back(0.f, 1.f, 0.f);
    vertices.emplace_back(0.5f, 0.5f, 0.f);
    vertices.emplace_back(0.51f, 0.5f, 0.f);
    
    std::vector<jtk::vec3<uint32_t>> triangles;
    triangles.emplace_back(0, 1, 4);
    triangles.emplace_back(1, 2, 5);
    triangles.emplace_back(2, 5, 3);
    triangles.emplace_back(3, 4, 0);
    
    TEST_ASSERT(!stitch_points(vertices, triangles, 0.0001f));
    TEST_ASSERT(stitch_points(vertices, triangles, 0.1f));
    TEST_EQ(5, vertices.size());
    }
  
  void test_undo_dyadic_subdivide()
    {
    std::vector<jtk::vec3<float>> vertices;
    vertices.emplace_back(0.f, 0.f, 0.f);
    vertices.emplace_back(1.f, 0.f, 0.f);
    vertices.emplace_back(1.f, 1.f, 0.f);
    vertices.emplace_back(0.f, 1.f, 0.f);
    vertices.emplace_back(0.5f, 0.5f, 0.f);
    std::vector<jtk::vec3<uint32_t>> triangles;
    triangles.emplace_back(0, 1, 4);
    triangles.emplace_back(1, 2, 4);
    triangles.emplace_back(2, 4, 3);
    triangles.emplace_back(3, 4, 0);

    jtk::dyadic_subdivide(vertices, triangles);
    jtk::dyadic_subdivide(vertices, triangles);
    TEST_EQ(64, triangles.size());

    jtk::undo_dyadic_subdivide(vertices, triangles);
    jtk::undo_dyadic_subdivide(vertices, triangles);
    TEST_EQ(4, triangles.size());
    TEST_EQ(5, vertices.size());
    }
  }

void run_all_geometry_tests()
  {
  using namespace jtk;
  test_stitch();
  test_undo_dyadic_subdivide();
  }
