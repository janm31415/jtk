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

  void test_delete_items()
    {
    std::vector<double> v;
    for (int i = 1; i <= 100; ++i)
      v.push_back((double)i);
    std::vector<int> indices_to_delete;
    indices_to_delete.push_back(9);
    indices_to_delete.push_back(19);
    indices_to_delete.push_back(29);
    indices_to_delete.push_back(39);
    indices_to_delete.push_back(49);
    indices_to_delete.push_back(59);
    indices_to_delete.push_back(69);
    indices_to_delete.push_back(79);
    indices_to_delete.push_back(89);
    indices_to_delete.push_back(99);
    delete_items(v, indices_to_delete);
    TEST_EQ(90, v.size());
    for (int i = 0; i < 90; ++i)
      {
      int value = (i / 9) * 10 + (i % 9) + 1;
      TEST_EQ((double)value, v[i]);
      }
    }
  }

void run_all_geometry_tests()
  {
  using namespace jtk;
  test_stitch();
  test_undo_dyadic_subdivide();
  test_delete_items();
  }
