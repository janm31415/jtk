#include "halffloat_tests.h"
#define JTK_HALFFLOAT_IMPLEMENTATION
#include "../jtk/halffloat.h"
#include "test_assert.h"


void run_all_halffloat_tests()
  {
  using namespace jtk;
  TEST_EQ(0.000000059604645f, halffloat_to_float(0b0000000000000001));

  TEST_EQ(0.000060975552f, halffloat_to_float(0b0000001111111111));

  TEST_EQ(0.000061035156f, halffloat_to_float(0b0000010000000000));

  TEST_EQ(65504.f, halffloat_to_float(0b0111101111111111));

  TEST_EQ(0.99951172f, halffloat_to_float(0b0011101111111111));

  TEST_EQ(1.f, halffloat_to_float(0b0011110000000000));

  TEST_EQ(1.00097656f, halffloat_to_float(0b0011110000000001));

  TEST_EQ(0.33325195f, halffloat_to_float(0b0011010101010101));

  TEST_EQ(-2.f, halffloat_to_float(0b1100000000000000));

  TEST_EQ(0.f, halffloat_to_float(0b0000000000000000));

  TEST_EQ(-0.f, halffloat_to_float(0b1000000000000000));

  TEST_EQ(0b0000000000000001, float_to_halffloat(halffloat_to_float(0b0000000000000001)));

  TEST_EQ(0b0000001111111111, float_to_halffloat(halffloat_to_float(0b0000001111111111)));

  TEST_EQ(0b0000010000000000, float_to_halffloat(halffloat_to_float(0b0000010000000000)));

  TEST_EQ(0b0111101111111111, float_to_halffloat(halffloat_to_float(0b0111101111111111)));

  TEST_EQ(0b0011101111111111, float_to_halffloat(halffloat_to_float(0b0011101111111111)));

  TEST_EQ(0b0011110000000000, float_to_halffloat(halffloat_to_float(0b0011110000000000)));

  TEST_EQ(0b0011110000000001, float_to_halffloat(halffloat_to_float(0b0011110000000001)));

  TEST_EQ(0b0011010101010101, float_to_halffloat(halffloat_to_float(0b0011010101010101)));

  TEST_EQ(0b1100000000000000, float_to_halffloat(halffloat_to_float(0b1100000000000000)));

  TEST_EQ(0b0000000000000000, float_to_halffloat(halffloat_to_float(0b0000000000000000)));

  TEST_EQ(0b1000000000000000, float_to_halffloat(halffloat_to_float(0b1000000000000000)));
  }
