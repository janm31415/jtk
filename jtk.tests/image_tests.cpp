#include "image_tests.h"

#include <stdint.h>

#define JTK_IMAGE_IMPLEMENTATION
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

#if !defined(_JTK_FOR_ARM)&&!defined(_JTK_NO_SIMD)
    inline image<uint32_t> census_transform_avx(const image<uint8_t>& im)
      {
      /*
      The census transform is an image operator that associates to each pixel of a grayscale image a binary string,
      encoding whether the pixel has smaller intensity than each of its neighbours, one for each bit.
      We use a 9x7 window.
      */
      using namespace image_details;
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
      using namespace image_details;
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
    
#endif
    }

  void test_census()
    {
#if !defined(_JTK_FOR_ARM)&&!defined(_JTK_NO_SIMD)
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
#endif
    }


  void test_census_32()
    {
#if !defined(_JTK_FOR_ARM)&&!defined(_JTK_NO_SIMD)
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
#endif
    }

  void test_image_integral()
    {
#if !defined(_JTK_FOR_ARM)&&!defined(_JTK_NO_SIMD)
    image<uint8_t> im(5, 5);
    im(0, 0) = 5;
    im(1, 0) = 2;
    im(2, 0) = 3;
    im(3, 0) = 4;
    im(4, 0) = 1;

    im(0, 1) = 1;
    im(1, 1) = 5;
    im(2, 1) = 4;
    im(3, 1) = 2;
    im(4, 1) = 3;

    im(0, 2) = 2;
    im(1, 2) = 2;
    im(2, 2) = 1;
    im(3, 2) = 3;
    im(4, 2) = 4;

    im(0, 3) = 3;
    im(1, 3) = 5;
    im(2, 3) = 6;
    im(3, 3) = 4;
    im(4, 3) = 5;

    im(0, 4) = 4;
    im(1, 4) = 1;
    im(2, 4) = 3;
    im(3, 4) = 2;
    im(4, 4) = 6;

    image<uint32_t> in;

    integral(in, im);

    TEST_EQ(6, in.width());
    TEST_EQ(6, in.height());

    TEST_EQ(0, in(0, 0));
    TEST_EQ(0, in(1, 0));
    TEST_EQ(0, in(2, 0));
    TEST_EQ(0, in(3, 0));
    TEST_EQ(0, in(4, 0));
    TEST_EQ(0, in(5, 0));

    TEST_EQ(0, in(0, 1));
    TEST_EQ(5, in(1, 1));
    TEST_EQ(7, in(2, 1));
    TEST_EQ(10, in(3, 1));
    TEST_EQ(14, in(4, 1));
    TEST_EQ(15, in(5, 1));

    TEST_EQ(0, in(0, 2));
    TEST_EQ(6, in(1, 2));
    TEST_EQ(13, in(2, 2));
    TEST_EQ(20, in(3, 2));
    TEST_EQ(26, in(4, 2));
    TEST_EQ(30, in(5, 2));

    TEST_EQ(0, in(0, 3));
    TEST_EQ(8, in(1, 3));
    TEST_EQ(17, in(2, 3));
    TEST_EQ(25, in(3, 3));
    TEST_EQ(34, in(4, 3));
    TEST_EQ(42, in(5, 3));

    TEST_EQ(0, in(0, 4));
    TEST_EQ(11, in(1, 4));
    TEST_EQ(25, in(2, 4));
    TEST_EQ(39, in(3, 4));
    TEST_EQ(52, in(4, 4));
    TEST_EQ(65, in(5, 4));

    TEST_EQ(0, in(0, 5));
    TEST_EQ(15, in(1, 5));
    TEST_EQ(30, in(2, 5));
    TEST_EQ(47, in(3, 5));
    TEST_EQ(62, in(4, 5));
    TEST_EQ(81, in(5, 5));
#endif
    }

  }


void run_all_image_tests()
  {
  using namespace jtk;
  test_fill_image();
  test_census();
  test_census_32();
  test_image_integral();
  }
