#include "image_tests.h"

#include <stdint.h>

#include "../jtk/image.h"
#include "../jtk/vec.h"
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
        TEST_EQ(5, im(x, y));
        }
      }
    fill_image(im, (uint8_t)6);
    for (int y = 0; y < h; ++y)
      {
      for (int x = 0; x < w; ++x)
        {
        TEST_EQ(6, im(x, y));
        }
      }
    }

  namespace
    {


    inline image<uint32_t> census_transform_avx(const image<uint8_t>& im)
      {
      /*
      The census transform is an image operator that associates to each pixel of a grayscale image a binary string,
      encoding whether the pixel has smaller intensity than each of its neighbours, one for each bit.
      We use a 9x7 window.
      */
      using namespace details;
      const int w = (int)im.width();
      const int s = (int)im.stride();
      const int h = (int)im.height();
      const int ns = (int)(sizeof(census_samples) / sizeof(int)) / 2;

      image<uint32_t> ct(w, h, true);

      const int cs = (int)ct.stride();

      for (int y = 4; y < h - 4; ++y)
        {
        __m256i descriptor = _mm256_setzero_si256();
        const uint8_t* p_im = im.data() + s * y;
        for (int x = 4; x < w - 4; ++x, ++p_im)
          {
          uint32_t* ct_pos = ct.data() + (y*cs) + x;
          uint8_t* descriptor_it = (uint8_t*)&descriptor;
          uint8_t value = *p_im;
          const __m256i value256 = _mm256_set1_epi8(value);
          const int pos = y * s + x;
          for (int p = 0; p < ns; p++)
            {
            *descriptor_it = *(p_im + census_samples[2 * p] + census_samples[2 * p + 1] * s);
            ++descriptor_it;
            }
          const __m256i r1 = _mm256_cmpgt_epi8(descriptor, value256);
          *ct_pos = _mm256_movemask_epi8(r1);
          }
        }
      return ct;
      }


    inline image<vec3<uint32_t>> census_transform_avx(const image<uint32_t>& im)
      {
      /*
      The census transform is an image operator that associates to each pixel of a grayscale image a binary string,
      encoding whether the pixel has smaller intensity than each of its neighbours, one for each bit.
      We use a 9x7 window.
      */
      using namespace details;
      const int w = (int)im.width();
      const int s = (int)im.stride();
      const int h = (int)im.height();
      const int ns = (int)(sizeof(census_samples) / sizeof(int)) / 2;

      image<vec3<uint32_t>> ct(w, h, true);

      const int cs = (int)ct.stride();

      for (int y = 4; y < h - 4; ++y)
        {
        __m256i descriptor_red = _mm256_setzero_si256();
        __m256i descriptor_green = _mm256_setzero_si256();
        __m256i descriptor_blue = _mm256_setzero_si256();
        const uint32_t* p_im = im.data() + s * y;
        for (int x = 4; x < w - 4; ++x, ++p_im)
          {
          vec3<uint32_t>* ct_pos = ct.data() + (y*cs) + x;
          uint8_t* descriptor_red_it = (uint8_t*)&descriptor_red;
          uint8_t* descriptor_green_it = (uint8_t*)&descriptor_green;
          uint8_t* descriptor_blue_it = (uint8_t*)&descriptor_blue;
          uint32_t value = *p_im;
          const __m256i red256 = _mm256_set1_epi8(value & 0xff);
          const __m256i green256 = _mm256_set1_epi8((value >> 8) & 0xff);
          const __m256i blue256 = _mm256_set1_epi8((value >> 16) & 0xff);
          const int pos = y * s + x;
          for (int p = 0; p < ns; p++)
            {
            const uint32_t color = *(p_im + census_samples[2 * p] + census_samples[2 * p + 1] * s);
            *descriptor_red_it = color & 0xff;
            *descriptor_green_it = (color >> 8) & 0xff;
            *descriptor_blue_it = (color >> 16) & 0xff;
            ++descriptor_red_it;
            ++descriptor_green_it;
            ++descriptor_blue_it;
            }
          const __m256i r1 = _mm256_cmpgt_epi8(descriptor_red, red256);
          const __m256i g1 = _mm256_cmpgt_epi8(descriptor_green, green256);
          const __m256i b1 = _mm256_cmpgt_epi8(descriptor_blue, blue256);
          *ct_pos = vec3<uint32_t>(_mm256_movemask_epi8(r1), _mm256_movemask_epi8(g1), _mm256_movemask_epi8(b1));
          }
        }
      return ct;
      }
    }

  void test_census()
    {
    int w = 64;
    int h = 64;
    image<uint8_t> im(w, h);
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        im(x, y) = (x*y) & 255;
    image<uint32_t> imc1 = census_transform_avx(im);
    image<uint32_t> imc2 = census_transform(im);
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        {
        TEST_EQ(imc1(x, y), imc2(x, y));
        }
    }


  void test_census_32()
    {
    int w = 64;
    int h = 64;
    image<uint32_t> im(w, h);
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        {
        uint32_t clr = ((x*y) & 255) | ((x*y*2) & 255) << 8 | ((x*y * 3) & 255) << 16 | 0xff000000;
        im(x, y) = clr;
        }
    image<vec3<uint32_t>> imc1 = census_transform_avx(im);
    image<uint32_t> imcr, imcg, imcb;
    census_transform(imcr, imcg, imcb, im);
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        {
        TEST_EQ(imcr(x, y), imc1(x, y)[0]);
        TEST_EQ(imcg(x, y), imc1(x, y)[1]);
        TEST_EQ(imcb(x, y), imc1(x, y)[2]);
        }
    }

  }


void run_all_image_tests()
  {
  using namespace jtk;
  test_fill_image();
  test_census();
  test_census_32();
  }
