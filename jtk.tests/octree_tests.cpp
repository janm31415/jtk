#include "octree_tests.h"

#include "../jtk/octree.h"

#include "test_assert.h"

namespace
  {

  uint32_t* crd(uint32_t a, uint32_t b, uint32_t c)
    {
    static uint32_t coord[3] = { 0, 0, 0 };
    coord[0] = a;
    coord[1] = b;
    coord[2] = c;
    return coord;
    }

  void test_indexed_octree_1()
    {
    jtk::indexed_octree<int32_t> oct(6);
    TEST_ASSERT(oct.tree_is_consistent());

    TEST_ASSERT(oct.find_leaf(crd(1, 2, 3)) == nullptr);
    uint8_t* child = oct.add_leaf(crd(1, 2, 3));
    TEST_ASSERT(oct.find_leaf(crd(1, 2, 3)) == child);
    oct.set_value(child, 3);
    TEST_EQ(3, oct.get_value(child));
    TEST_ASSERT(oct.tree_is_consistent());

    uint8_t* child2 = oct.add_leaf(crd(63, 63, 63));
    TEST_ASSERT(oct.find_leaf(crd(63, 63, 63)) == child2);
    TEST_ASSERT(child != child2);
    TEST_ASSERT(oct.tree_is_consistent());

    TEST_ASSERT(oct.in_bounds(crd(63, 63, 63), oct.max_depth()));
    TEST_ASSERT(!oct.in_bounds(crd(64, 63, 63), oct.max_depth()));

    oct.for_each_leaf([&](uint8_t* node, uint8_t* /*parent*/, uint32_t* /*coord*/)
      {
      oct.set_value(node, 5);
      });

    TEST_EQ(5, oct.get_value(child));
    TEST_EQ(5, oct.get_value(child2));

    oct.remove_leaf(crd(63, 63, 63));
    TEST_ASSERT(oct.tree_is_consistent());
    TEST_ASSERT(oct.find_leaf(crd(63, 63, 63)) == nullptr);
    TEST_ASSERT(oct.find_leaf(crd(1, 2, 3)) == child);
    }

  void test_octree_1()
    {
    jtk::octree<int32_t> oct(6);
    TEST_ASSERT(oct.tree_is_consistent());

    jtk::vec3<float> pt1(0.1f,0.02f,0.3f);
    jtk::vec3<float> pt2(0.99f, 0.99f, 0.99f);

    TEST_ASSERT(oct.find_leaf(pt1) == nullptr);
    uint8_t* child = oct.add_leaf(pt1);
    TEST_ASSERT(oct.find_leaf(pt1) == child);
    oct.set_value(child, 3);
    TEST_EQ(3, oct.get_value(child));
    TEST_ASSERT(oct.tree_is_consistent());
    
    uint8_t* child2 = oct.add_leaf(pt2);
    TEST_ASSERT(oct.find_leaf(pt2) == child2);
    TEST_ASSERT(child != child2);
    TEST_ASSERT(oct.tree_is_consistent());

    oct.for_each_leaf([&](uint8_t* leaf, const jtk::vec3<double>& center, double halfsize)
      {
      TEST_EQ(0.5/(double)(1<<oct.max_depth()), halfsize);
      oct.set_value(leaf, 107);
      if (leaf == child2)
        {
        TEST_EQ(1.0 - halfsize, center[0]);
        TEST_EQ(1.0 - halfsize, center[1]);
        TEST_EQ(1.0 - halfsize, center[2]);
        }
      });

    TEST_EQ(107, oct.get_value(child));
    TEST_EQ(107, oct.get_value(child2));

    oct.remove_leaf(pt2);
    TEST_ASSERT(oct.tree_is_consistent());
    TEST_ASSERT(oct.find_leaf(pt2) == nullptr);
    TEST_ASSERT(oct.find_leaf(pt1) == child);
    }
  }

void run_all_octree_tests()
  {
  test_indexed_octree_1();
  test_octree_1();
  }