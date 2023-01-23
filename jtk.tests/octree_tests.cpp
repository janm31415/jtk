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
    }
  }

void run_all_octree_tests()
  {
  test_indexed_octree_1();
  }