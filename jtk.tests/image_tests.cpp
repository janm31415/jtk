#include "image_tests.h"

#include <stdint.h>

#include "../jtk/image.h"
#include "test_assert.h"

namespace jtk
  {

  void test_fill_image()
    {
    int w = 10;
    int h = 10;
    image<uint8_t> im(w, h);
    fill_image(im, (uint8_t)5);
    for (int y = 0; y < h; ++y)
      {
      for (int x = 0; x < w; ++x)
        {
        TEST_EQ(5, im(x,y));
        }
      }
    fill_image(im, (uint8_t)6);
    for (int y = 0; y < h; ++y)
      {
      for (int x = 0; x < w; ++x)
        {
        TEST_EQ(6, im(x,y));
        }
      }
    }
  
  }


void run_all_image_tests()
  {
  using namespace jtk;
  test_fill_image();
  }
