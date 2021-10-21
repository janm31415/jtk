#include "qbvh_tests.h"

#define JTK_QBVH_IMPLEMENTATION
#include "../jtk/qbvh.h"
#include "test_assert.h"

#include <iostream>
#include <chrono>
#include <random>

namespace jtk
  {
  void test_float4_1()
    {
    float4 f;
    f[0] = 2.f;
    f[1] = 3.f;
    f[2] = 4.f;
    f[3] = 5.f;
    TEST_EQ(f[0], 2.f);
    TEST_EQ(f[1], 3.f);
    TEST_EQ(f[2], 4.f);
    TEST_EQ(f[3], 5.f);
    }

  void test_float4_cross()
    {
    float4 a(1.f, 0.f, 0.f, 1.f);
    float4 b(0.f, 1.f, 0.f, 1.f);
    float4 c = cross(a, b);
    TEST_EQ(0.f, c[0]);
    TEST_EQ(0.f, c[1]);
    TEST_EQ(1.f, c[2]);
    TEST_EQ(0.f, c[3]);
    }

  void test_float4_cross2()
    {
    float4 a(1.f, 2.f, 3.f, 1.f);
    float4 b(-1.f, 0.0f, -0.5f, 1.f);
    float4 c = cross(a, b);
    TEST_EQ(-1.f, c[0]);
    TEST_EQ(-2.5f, c[1]);
    TEST_EQ(2.f, c[2]);
    TEST_EQ(0.f, c[3]);
    }

  void test_float4_dot()
    {
    float4 a(1.f, 2.f, 3.f, 1.f);
    float4 b(-1.f, 0.0f, -0.5f, 1.f);
    float d = dot(a, b);
    TEST_EQ(-1.f - 0.5*3.f, d);
    }

  void test_float4_rsqrt()
    {
    float4 a(1.f, 2.f, 3.f, 4.f);
    auto res = rsqrt(a);
    float tol = 1e-7f;
    TEST_EQ_CLOSE(1.f, res[0], tol);
    TEST_EQ_CLOSE(1.f / std::sqrt(2.f), res[1], tol);
    TEST_EQ_CLOSE(1.f / std::sqrt(3.f), res[2], tol);
    TEST_EQ_CLOSE(0.5f, res[3], tol);

    float4 b(1.f, 0.f, 0.f, 4.f);
    res = rsqrt(b);
    TEST_EQ_CLOSE(1.f, res[0], tol);
    TEST_EQ_CLOSE(std::numeric_limits<float>::infinity(), res[1], tol);
    TEST_EQ_CLOSE(std::numeric_limits<float>::infinity(), res[2], tol);
    TEST_EQ_CLOSE(0.5f, res[3], tol);
    }

  void test_float4()
    {
    float4 one(1.f), two(2.f), iota(1.f, 2.f, 3.f, 4.f);
    TEST_EQ(1.f, min_horizontal(iota));
    TEST_EQ(4.f, max_horizontal(iota));
    TEST_EQ(one[0], 1.f);
    TEST_EQ(one[1], 1.f);
    TEST_EQ(one[2], 1.f);
    TEST_EQ(one[3], 1.f);
    TEST_EQ(two[0], 2.f);
    TEST_EQ(two[1], 2.f);
    TEST_EQ(two[2], 2.f);
    TEST_EQ(two[3], 2.f);
    TEST_EQ(iota[0], 1.f);
    TEST_EQ(iota[1], 2.f);
    TEST_EQ(iota[2], 3.f);
    TEST_EQ(iota[3], 4.f);
    float4 res1 = one + iota;
    TEST_EQ(res1[0], 2.f);
    TEST_EQ(res1[1], 3.f);
    TEST_EQ(res1[2], 4.f);
    TEST_EQ(res1[3], 5.f);
    float4 res2 = two - iota;
    TEST_EQ(res2[0], 1.f);
    TEST_EQ(res2[1], 0.f);
    TEST_EQ(res2[2], -1.f);
    TEST_EQ(res2[3], -2.f);
    TEST_EQ(-two[0], -2.f);
    TEST_EQ(-two[1], -2.f);
    TEST_EQ(-two[2], -2.f);
    TEST_EQ(-two[3], -2.f);
    float4 res3 = two * iota;
    TEST_EQ(res3[0], 2.f);
    TEST_EQ(res3[1], 4.f);
    TEST_EQ(res3[2], 6.f);
    TEST_EQ(res3[3], 8.f);
    float4 res4 = 5 * iota;
    TEST_EQ(res4[0], 5.f);
    TEST_EQ(res4[1], 10.f);
    TEST_EQ(res4[2], 15.f);
    TEST_EQ(res4[3], 20.f);
    float4 res5 = iota * 6;
    TEST_EQ(res5[0], 6.f);
    TEST_EQ(res5[1], 12.f);
    TEST_EQ(res5[2], 18.f);
    TEST_EQ(res5[3], 24.f);
    float4 res6 = min(two, iota);
    TEST_EQ(res6[0], 1.f);
    TEST_EQ(res6[1], 2.f);
    TEST_EQ(res6[2], 2.f);
    TEST_EQ(res6[3], 2.f);
    float4 res7 = max(two, iota);
    TEST_EQ(res7[0], 2.f);
    TEST_EQ(res7[1], 2.f);
    TEST_EQ(res7[2], 3.f);
    TEST_EQ(res7[3], 4.f);
    bool4 res8 = (two == iota);
    TEST_EQ(res8[0], 0);
    TEST_EQ(res8[1], -1);
    TEST_EQ(res8[2], 0);
    TEST_EQ(res8[3], 0);
    bool4 res9 = (two < iota);
    TEST_EQ(res9[0], 0);
    TEST_EQ(res9[1], 0);
    TEST_EQ(res9[2], -1);
    TEST_EQ(res9[3], -1);
    float4 res10 = reciprocal(float4(0.f, 1.f, 2.f, 3.f));
    float tol = 1e-7f;
    TEST_EQ(res10[0], std::numeric_limits<float>::infinity());
    TEST_EQ_CLOSE(res10[1], 1.f, tol);
    TEST_EQ_CLOSE(res10[2], 0.5f, tol);
    TEST_EQ_CLOSE(res10[3], 0.333333333333333f, tol);
    float4 res11 = masked_update(bool4(true, false, false, true), two, iota);
    TEST_EQ(res11[0], 1.f);
    TEST_EQ(res11[1], 2.f);
    TEST_EQ(res11[2], 2.f);
    TEST_EQ(res11[3], 4.f);
    }

  void test_float4_transpose()
    {
    float4 c0(1.f, 2.f, 3.f, 4.f);
    float4 c1(5.f, 6.f, 7.f, 8.f);
    float4 c2(9.f, 10.f, 11.f, 12.f);
    float4 c3(13.f, 14.f, 15.f, 16.f);
    float4 r0, r1, r2, r3;
    transpose(r0, r1, r2, r3, c0, c1, c2, c3);
    TEST_EQ(r0[0], 1.f);
    TEST_EQ(r0[1], 5.f);
    TEST_EQ(r0[2], 9.f);
    TEST_EQ(r0[3], 13.f);
    TEST_EQ(r1[0], 2.f);
    TEST_EQ(r1[1], 6.f);
    TEST_EQ(r1[2], 10.f);
    TEST_EQ(r1[3], 14.f);
    TEST_EQ(r2[0], 3.f);
    TEST_EQ(r2[1], 7.f);
    TEST_EQ(r2[2], 11.f);
    TEST_EQ(r2[3], 15.f);
    TEST_EQ(r3[0], 4.f);
    TEST_EQ(r3[1], 8.f);
    TEST_EQ(r3[2], 12.f);
    TEST_EQ(r3[3], 16.f);
    }

  void test_bool4()
    {
    bool4 b(true, false, false, true);
    bool4 b2(true);
    bool4 b3(false);
    TEST_EQ(b[0], -1);
    TEST_EQ(b[1], 0);
    TEST_EQ(b[2], 0);
    TEST_EQ(b[3], -1);
    TEST_ASSERT(!all(b));
    TEST_ASSERT(all(b2));
    TEST_ASSERT(!all(b3));
    TEST_ASSERT(!none(b2));
    TEST_ASSERT(none(b3));
    bool4 notb = !b;
    TEST_EQ(notb[0], 0);
    TEST_EQ(notb[1], -1);
    TEST_EQ(notb[2], -1);
    TEST_EQ(notb[3], 0);
    bool4 b4 = bool4(true, false, false, false) | bool4(false, true, true, false);
    TEST_EQ(b4[0], -1);
    TEST_EQ(b4[1], -1);
    TEST_EQ(b4[2], -1);
    TEST_EQ(b4[3], 0);
    bool4 b5 = bool4(true, false, false, false) & bool4(true, true, true, false);
    TEST_EQ(b5[0], -1);
    TEST_EQ(b5[1], 0);
    TEST_EQ(b5[2], 0);
    TEST_EQ(b5[3], 0);
    bool4 b6 = bool4(true, false, true, false) ^ bool4(false, true, true, false);
    TEST_EQ(b6[0], -1);
    TEST_EQ(b6[1], -1);
    TEST_EQ(b6[2], 0);
    TEST_EQ(b6[3], 0);
    TEST_ASSERT(none(b3&b2));
    TEST_ASSERT(all(b3 | b2));
    bool4 b7 = bool4(true, false, true, false) == bool4(false, true, true, false);
    TEST_EQ(b7[0], 0);
    TEST_EQ(b7[1], 0);
    TEST_EQ(b7[2], -1);
    TEST_EQ(b7[3], -1);
    bool4 b8 = bool4(true, false, true, false) != bool4(false, true, true, false);
    TEST_EQ(b8[0], -1);
    TEST_EQ(b8[1], -1);
    TEST_EQ(b8[2], -0);
    TEST_EQ(b8[3], -0);
    }

  void test_int4()
    {
    int4 one(1), two(2), iota(1, 2, 3, 4);
    TEST_EQ(one[0], 1);
    TEST_EQ(one[1], 1);
    TEST_EQ(one[2], 1);
    TEST_EQ(one[3], 1);
    TEST_EQ(two[0], 2);
    TEST_EQ(two[1], 2);
    TEST_EQ(two[2], 2);
    TEST_EQ(two[3], 2);
    TEST_EQ(iota[0], 1);
    TEST_EQ(iota[1], 2);
    TEST_EQ(iota[2], 3);
    TEST_EQ(iota[3], 4);
    int4 res1 = one + iota;
    TEST_EQ(res1[0], 2);
    TEST_EQ(res1[1], 3);
    TEST_EQ(res1[2], 4);
    TEST_EQ(res1[3], 5);
    int4 res2 = two - iota;
    TEST_EQ(res2[0], 1);
    TEST_EQ(res2[1], 0);
    TEST_EQ(res2[2], -1);
    TEST_EQ(res2[3], -2);
    TEST_EQ(-two[0], -2);
    TEST_EQ(-two[1], -2);
    TEST_EQ(-two[2], -2);
    TEST_EQ(-two[3], -2);
    int4 res3 = two * iota;
    TEST_EQ(res3[0], 2);
    TEST_EQ(res3[1], 4);
    TEST_EQ(res3[2], 6);
    TEST_EQ(res3[3], 8);
    int4 res4 = 5 * iota;
    TEST_EQ(res4[0], 5);
    TEST_EQ(res4[1], 10);
    TEST_EQ(res4[2], 15);
    TEST_EQ(res4[3], 20);
    int4 res5 = iota * 6;
    TEST_EQ(res5[0], 6);
    TEST_EQ(res5[1], 12);
    TEST_EQ(res5[2], 18);
    TEST_EQ(res5[3], 24);
    int4 res6 = min(two, iota);
    TEST_EQ(res6[0], 1);
    TEST_EQ(res6[1], 2);
    TEST_EQ(res6[2], 2);
    TEST_EQ(res6[3], 2);
    int4 res7 = max(two, iota);
    TEST_EQ(res7[0], 2);
    TEST_EQ(res7[1], 2);
    TEST_EQ(res7[2], 3);
    TEST_EQ(res7[3], 4);
    bool4 res8 = (two == iota);
    TEST_EQ(res8[0], 0);
    TEST_EQ(res8[1], -1);
    TEST_EQ(res8[2], 0);
    TEST_EQ(res8[3], 0);
    bool4 res9 = (two < iota);
    TEST_EQ(res9[0], 0);
    TEST_EQ(res9[1], 0);
    TEST_EQ(res9[2], -1);
    TEST_EQ(res9[3], -1);
#ifndef _JTK_FOR_ARM
    int4 res10 = (iota << 2);
    TEST_EQ(res10[0], 4);
    TEST_EQ(res10[1], 8);
    TEST_EQ(res10[2], 12);
    TEST_EQ(res10[3], 16);
#endif
    int4 res11 = masked_update(bool4(true, false, false, true), two, iota);
    TEST_EQ(res11[0], 1);
    TEST_EQ(res11[1], 2);
    TEST_EQ(res11[2], 2);
    TEST_EQ(res11[3], 4);
    }


  namespace
    {

    inline void invert(float* out, const float* m)
      {
      float inv[16], det;
      int i;

      inv[0] = m[5] * m[10] * m[15] -
        m[5] * m[11] * m[14] -
        m[9] * m[6] * m[15] +
        m[9] * m[7] * m[14] +
        m[13] * m[6] * m[11] -
        m[13] * m[7] * m[10];

      inv[4] = -m[4] * m[10] * m[15] +
        m[4] * m[11] * m[14] +
        m[8] * m[6] * m[15] -
        m[8] * m[7] * m[14] -
        m[12] * m[6] * m[11] +
        m[12] * m[7] * m[10];

      inv[8] = m[4] * m[9] * m[15] -
        m[4] * m[11] * m[13] -
        m[8] * m[5] * m[15] +
        m[8] * m[7] * m[13] +
        m[12] * m[5] * m[11] -
        m[12] * m[7] * m[9];

      inv[12] = -m[4] * m[9] * m[14] +
        m[4] * m[10] * m[13] +
        m[8] * m[5] * m[14] -
        m[8] * m[6] * m[13] -
        m[12] * m[5] * m[10] +
        m[12] * m[6] * m[9];

      inv[1] = -m[1] * m[10] * m[15] +
        m[1] * m[11] * m[14] +
        m[9] * m[2] * m[15] -
        m[9] * m[3] * m[14] -
        m[13] * m[2] * m[11] +
        m[13] * m[3] * m[10];

      inv[5] = m[0] * m[10] * m[15] -
        m[0] * m[11] * m[14] -
        m[8] * m[2] * m[15] +
        m[8] * m[3] * m[14] +
        m[12] * m[2] * m[11] -
        m[12] * m[3] * m[10];

      inv[9] = -m[0] * m[9] * m[15] +
        m[0] * m[11] * m[13] +
        m[8] * m[1] * m[15] -
        m[8] * m[3] * m[13] -
        m[12] * m[1] * m[11] +
        m[12] * m[3] * m[9];

      inv[13] = m[0] * m[9] * m[14] -
        m[0] * m[10] * m[13] -
        m[8] * m[1] * m[14] +
        m[8] * m[2] * m[13] +
        m[12] * m[1] * m[10] -
        m[12] * m[2] * m[9];

      inv[2] = m[1] * m[6] * m[15] -
        m[1] * m[7] * m[14] -
        m[5] * m[2] * m[15] +
        m[5] * m[3] * m[14] +
        m[13] * m[2] * m[7] -
        m[13] * m[3] * m[6];

      inv[6] = -m[0] * m[6] * m[15] +
        m[0] * m[7] * m[14] +
        m[4] * m[2] * m[15] -
        m[4] * m[3] * m[14] -
        m[12] * m[2] * m[7] +
        m[12] * m[3] * m[6];

      inv[10] = m[0] * m[5] * m[15] -
        m[0] * m[7] * m[13] -
        m[4] * m[1] * m[15] +
        m[4] * m[3] * m[13] +
        m[12] * m[1] * m[7] -
        m[12] * m[3] * m[5];

      inv[14] = -m[0] * m[5] * m[14] +
        m[0] * m[6] * m[13] +
        m[4] * m[1] * m[14] -
        m[4] * m[2] * m[13] -
        m[12] * m[1] * m[6] +
        m[12] * m[2] * m[5];

      inv[3] = -m[1] * m[6] * m[11] +
        m[1] * m[7] * m[10] +
        m[5] * m[2] * m[11] -
        m[5] * m[3] * m[10] -
        m[9] * m[2] * m[7] +
        m[9] * m[3] * m[6];

      inv[7] = m[0] * m[6] * m[11] -
        m[0] * m[7] * m[10] -
        m[4] * m[2] * m[11] +
        m[4] * m[3] * m[10] +
        m[8] * m[2] * m[7] -
        m[8] * m[3] * m[6];

      inv[11] = -m[0] * m[5] * m[11] +
        m[0] * m[7] * m[9] +
        m[4] * m[1] * m[11] -
        m[4] * m[3] * m[9] -
        m[8] * m[1] * m[7] +
        m[8] * m[3] * m[5];

      inv[15] = m[0] * m[5] * m[10] -
        m[0] * m[6] * m[9] -
        m[4] * m[1] * m[10] +
        m[4] * m[2] * m[9] +
        m[8] * m[1] * m[6] -
        m[8] * m[2] * m[5];

      det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

      for (i = 0; i < 16; ++i)
        out[i] = inv[i] / det;
      }

    inline void invert_orthonormal(float* out, const float* in)
      {
      // transpose
      out[0] = in[0];
      out[1] = in[4];
      out[2] = in[8];
      out[4] = in[1];
      out[5] = in[5];
      out[6] = in[9];
      out[8] = in[2];
      out[9] = in[6];
      out[10] = in[10];
      out[3] = 0;
      out[7] = 0;
      out[11] = 0;
      out[15] = 1;

      // translation
      out[12] = -(in[0] * in[12] + in[1] * in[13] + in[2] * in[14]);
      out[13] = -(in[4] * in[12] + in[5] * in[13] + in[6] * in[14]);
      out[14] = -(in[8] * in[12] + in[9] * in[13] + in[10] * in[14]);
      }

    inline void matrix_multiply(float* out, const float* left, const float* right)
      {
      for (int i = 0; i < 4; ++i)
        {
        for (int j = 0; j < 4; ++j)
          {
          out[i + (j << 2)] = left[i] * right[(j << 2)] + left[i + 4] * right[(j << 2) + 1] + left[i + 8] * right[(j << 2) + 2] + left[i + 12] * right[(j << 2) + 3];
          }
        }
      }
    }

  void test_float4x4()
    {
    float4x4 m, m2, m_inv, m_inv_inv;
    float inv[16];
    m = get_identity();

    TEST_EQ(1.f, m[0]);
    TEST_EQ(0.f, m[1]);
    TEST_EQ(0.f, m[2]);
    TEST_EQ(0.f, m[3]);
    TEST_EQ(0.f, m[4]);
    TEST_EQ(1.f, m[5]);
    TEST_EQ(0.f, m[6]);
    TEST_EQ(0.f, m[7]);
    TEST_EQ(0.f, m[8]);
    TEST_EQ(0.f, m[9]);
    TEST_EQ(1.f, m[10]);
    TEST_EQ(0.f, m[11]);
    TEST_EQ(0.f, m[12]);
    TEST_EQ(0.f, m[13]);
    TEST_EQ(0.f, m[14]);
    TEST_EQ(1.f, m[15]);

    for (int i = 0; i < 16; ++i)
      m2[i] = (float)i;
    m = transpose(m2);
    TEST_EQ(0.f, m[0]);
    TEST_EQ(4.f, m[1]);
    TEST_EQ(8.f, m[2]);
    TEST_EQ(12.f, m[3]);
    TEST_EQ(1.f, m[4]);
    TEST_EQ(5.f, m[5]);
    TEST_EQ(9.f, m[6]);
    TEST_EQ(13.f, m[7]);
    TEST_EQ(2.f, m[8]);
    TEST_EQ(6.f, m[9]);
    TEST_EQ(10.f, m[10]);
    TEST_EQ(14.f, m[11]);
    TEST_EQ(3.f, m[12]);
    TEST_EQ(7.f, m[13]);
    TEST_EQ(11.f, m[14]);
    TEST_EQ(15.f, m[15]);

    m = get_identity();
    m[0] = 1.f / std::sqrt(2.f);
    m[1] = 1.f / std::sqrt(2.f);
    m[4] = -1.f / std::sqrt(2.f);
    m[5] = 1.f / std::sqrt(2.f);
    m[10] = 1.f;
    m[12] = 31.f;
    m[13] = 52.f;
    m[14] = -27.f;

    invert_orthonormal(inv, m.f);

    m_inv = invert_orthonormal(m);

    for (int i = 0; i < 16; ++i)
      {
      TEST_EQ(inv[i], m_inv[i]);
      }

    m_inv_inv = invert_orthonormal(m_inv);
    for (int i = 0; i < 16; ++i)
      {
      TEST_EQ_CLOSE(m[i], m_inv_inv[i], 1e-5f);
      }

    m = get_identity();
    m[0] = 1.f / std::sqrt(2.f);
    m[1] = 1.f / std::sqrt(2.f);
    m[4] = -1.f / std::sqrt(2.f);
    m[5] = 1.f / std::sqrt(2.f);
    m[10] = 1.f;
    m[12] = 31.f;
    m[13] = 52.f;
    m[14] = -27.f;

    invert(inv, m.f);

    m_inv = invert(m);

    for (int i = 0; i < 16; ++i)
      {
      TEST_EQ_CLOSE(inv[i], m_inv[i], 1e-5f);
      }

    m_inv_inv = invert(m_inv);
    for (int i = 0; i < 16; ++i)
      {
      TEST_EQ_CLOSE(m[i], m_inv_inv[i], 1e-5f);
      }

    m[0] = 0.f;
    m[1] = 3.f;
    m[2] = 0.f;
    m[3] = 0.f;

    m[4] = 2.f;
    m[5] = 7.f;
    m[6] = 0.f;
    m[7] = 1.f;

    m[8] = 8.f;
    m[9] = 1.f;
    m[10] = 1.f;
    m[11] = 0.f;

    m[12] = 6.f;
    m[13] = 0.f;
    m[14] = 2.f;
    m[15] = 1.f;

    invert(inv, m.f);

    m_inv = invert(m);

    for (int i = 0; i < 16; ++i)
      {
      TEST_EQ_CLOSE(inv[i], m_inv[i], 1e-5f);
      }

    m_inv_inv = invert(m_inv);
    for (int i = 0; i < 16; ++i)
      {
      TEST_EQ_CLOSE(m[i], m_inv_inv[i], 1e-5f);
      }

    m = get_identity();
    m[0] = 1.f / std::sqrt(2.f);
    m[1] = 1.f / std::sqrt(2.f);
    m[4] = -1.f / std::sqrt(2.f);
    m[5] = 1.f / std::sqrt(2.f);
    m[10] = 1.f;
    m[12] = 31.f;
    m[13] = 52.f;
    m[14] = -27.f;

    float4 v0(1.f, 0.f, 0.f, 0.f);
    auto res = matrix_vector_multiply(m, v0);
    TEST_EQ(m[0], res[0]);
    TEST_EQ(m[1], res[1]);
    TEST_EQ(m[2], res[2]);
    TEST_EQ(m[3], res[3]);

    float4 v1(0.f, 1.f, 0.f, 0.f);
    res = matrix_vector_multiply(m, v1);
    TEST_EQ(m[4], res[0]);
    TEST_EQ(m[5], res[1]);
    TEST_EQ(m[6], res[2]);
    TEST_EQ(m[7], res[3]);

    float4 v2(0.f, 0.f, 1.f, 0.f);
    res = matrix_vector_multiply(m, v2);
    TEST_EQ(m[8], res[0]);
    TEST_EQ(m[9], res[1]);
    TEST_EQ(m[10], res[2]);
    TEST_EQ(m[11], res[3]);

    float4 v3(0.f, 0.f, 0.f, 1.f);
    res = matrix_vector_multiply(m, v3);
    TEST_EQ(m[12], res[0]);
    TEST_EQ(m[13], res[1]);
    TEST_EQ(m[14], res[2]);
    TEST_EQ(m[15], res[3]);

    m = get_identity();
    m[0] = 1.f / std::sqrt(2.f);
    m[1] = 1.f / std::sqrt(2.f);
    m[4] = -1.f / std::sqrt(2.f);
    m[5] = 1.f / std::sqrt(2.f);
    m[10] = 1.f;
    m[12] = 31.f;
    m[13] = 52.f;
    m[14] = -27.f;

    m_inv = invert_orthonormal(m);
    m2 = matrix_matrix_multiply(m, m_inv);

    TEST_EQ_CLOSE(1.f, m2[0], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[1], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[2], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[3], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[4], 1e-5f);
    TEST_EQ_CLOSE(1.f, m2[5], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[6], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[7], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[8], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[9], 1e-5f);
    TEST_EQ_CLOSE(1.f, m2[10], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[11], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[12], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[13], 1e-5f);
    TEST_EQ_CLOSE(0.f, m2[14], 1e-5f);
    TEST_EQ_CLOSE(1.f, m2[15], 1e-5f);
    }

  void run_all_simd_tests()
    {
    test_float4_1();
    test_float4_cross();
    test_float4_cross2();
    test_float4_dot();
    test_float4_rsqrt();
    test_float4();
    test_float4_transpose();
    test_bool4();
    test_int4();
    test_float4x4();
    }

  struct fixture_cube
    {
    std::vector<vec3<float>> vertices;
    std::vector<vec3<uint32_t>> triangles;

    fixture_cube()
      {
      vertices.emplace_back(-1.f, -1.f, 1.f);
      vertices.emplace_back(1.f, -1.f, 1.f);
      vertices.emplace_back(1.f, 1.f, 1.f);
      vertices.emplace_back(-1.f, 1.f, 1.f);
      vertices.emplace_back(-1.f, -1.f, -1.f);
      vertices.emplace_back(1.f, -1.f, -1.f);
      vertices.emplace_back(1.f, 1.f, -1.f);
      vertices.emplace_back(-1.f, 1.f, -1.f);
      triangles.emplace_back(0, 1, 2);
      triangles.emplace_back(0, 2, 3);
      triangles.emplace_back(7, 6, 5);
      triangles.emplace_back(7, 5, 4);
      triangles.emplace_back(1, 0, 4);
      triangles.emplace_back(1, 4, 5);
      triangles.emplace_back(2, 1, 5);
      triangles.emplace_back(2, 5, 6);
      triangles.emplace_back(3, 2, 6);
      triangles.emplace_back(3, 6, 7);
      triangles.emplace_back(0, 3, 7);
      triangles.emplace_back(0, 7, 4);
      }
    };

  struct find_closest_with_ray : public fixture_cube {
    void test()
      {
      uint32_t nr_of_triangles = (uint32_t)triangles.size();
      qbvh_voxel total_bb, centroid_bb;
      auto qbvh_voxels = build_triangle_qbvh_voxels(total_bb, centroid_bb, vertices.data(), triangles.data(), nr_of_triangles);
      qbvh::properties props;
      props.leaf_size = 2;
      qbvh bvh(qbvh_voxels, nr_of_triangles, total_bb, centroid_bb, props);
      delete[] qbvh_voxels;

      float4 origin(0.8f, -0.5f, 0.8f, 1.f);
      float4 dir(1.0f, 0.0f, 0.f, 0.f);

      ray r;
      r.orig = origin;
      r.dir = dir;
      r.t_near = 0.f;
      r.t_far = std::numeric_limits<float>::max();

      uint32_t triangle_id;

      auto hit = bvh.find_closest_triangle(triangle_id, r, triangles.data(), vertices.data());

      TEST_ASSERT(std::abs(hit.distance - 0.2f) <= 1e-6f);
      TEST_ASSERT(triangle_id == 6); // v1 v2 v5

      r.t_near = -std::numeric_limits<float>::max();
      r.t_far = 0.f;
      hit = bvh.find_closest_triangle(triangle_id, r, triangles.data(), vertices.data());
      TEST_ASSERT(std::abs(hit.distance + 1.8f) <= 1e-6f);
      TEST_ASSERT(triangle_id == 10); // v0 v3 v7

      r.t_near = -std::numeric_limits<float>::max();
      r.t_far = std::numeric_limits<float>::max();
      hit = bvh.find_closest_triangle(triangle_id, r, triangles.data(), vertices.data());

      TEST_ASSERT(std::abs(hit.distance - 0.2f) <= 1e-6f);
      TEST_ASSERT(triangle_id == 6); // v1 v2 v5

      }
    };

  void run_all_bvh_tests()
    {
    find_closest_with_ray().test();
    }

#ifdef _DEBUG
#define VEC_SIZE 1000
#else
#define VEC_SIZE 100000
#endif

  void test_std_partition()
    {
    std::mt19937 gen(0);
    std::uniform_int_distribution<> dis(-10, 10);
    std::vector<int> number(VEC_SIZE);
    for (auto& n : number)
      n = dis(gen);
    std::vector<int> number_copy = number;
    uint64_t total_time = 0;
    std::vector<int>::iterator mid;
    for (int iter = 0; iter < 1000; ++iter)
      {
      number = number_copy;
      auto tic = std::chrono::high_resolution_clock::now();
      mid = std::partition(number.begin(), number.end(), [](int n) { return n < 0; });
      auto toc = std::chrono::high_resolution_clock::now();
      if (iter > 500)
        total_time += std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count();
      }
    for (auto it = number.begin(); it != mid; ++it)
      TEST_ASSERT(*it < 0);
    for (auto it = mid; it != number.end(); ++it)
      TEST_ASSERT(*it >= 0);
    std::cout << "std_partition: " << total_time << "us" << std::endl;
    }

  void test_parallel_partition()
    {
    std::mt19937 gen(0);
    std::uniform_int_distribution<> dis(-10, 10);
    std::vector<int> number(VEC_SIZE);
    for (auto& n : number)
      n = dis(gen);
    std::vector<int> number_copy = number;
    uint64_t total_time = 0;
    std::vector<int>::iterator mid;
    for (int iter = 0; iter < 1000; ++iter)
      {
      number = number_copy;
      auto tic = std::chrono::high_resolution_clock::now();
      mid = parallel_partition(number.begin(), number.end(), [](int n) { return n < 0; });
      auto toc = std::chrono::high_resolution_clock::now();
      if (iter > 500)
        total_time += std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count();
      }
    for (auto it = number.begin(); it != mid; ++it)
      TEST_ASSERT(*it < 0);
    for (auto it = mid; it != number.end(); ++it)
      TEST_ASSERT(*it >= 0);
    std::cout << "parallel_partition: " << total_time << "us" << std::endl;
    }

  void run_all_parallel_partition_tests()
    {
    test_std_partition();
    test_parallel_partition();
    }

  void aligned_vector_test_1()
    {
    aligned_vector<double> v(5);
    TEST_EQ(5, v.size());
    v.push_back(3.0);
    TEST_EQ(3.0, v[5]);
    aligned_vector<double> v2;
    TEST_ASSERT(v2.empty());
    }

  void aligned_vector_test_2()
    {
    aligned_vector<double> v;
    v.reserve(100);
    TEST_EQ(100, v.capacity());
    }

  void aligned_vector_test_3()
    {
    aligned_vector<double> v;
    v.push_back(1.0);
    v.push_back(2.0);
    v.push_back(3.0);
    v.push_back(4.0);
    v.push_back(5.0);
    auto it = v.begin();
    auto it_end = v.end();
    double val = 1.0;
    for (; it != it_end; ++it)
      {
      TEST_EQ(val, *it);
      val += 1.0;
      }

    aligned_vector<double> v2(v);
    it = v.begin();
    val = 1.0;
    for (; it != it_end; ++it)
      {
      TEST_EQ(val, *it);
      val += 1.0;
      }

    for (auto& value : v2)
      value *= 2.0;
    it = v2.begin();
    it_end = v2.end();
    val = 1.0;
    for (; it != it_end; ++it)
      {
      TEST_EQ(val*2.0, *it);
      val += 1.0;
      }
    }

  void aligned_vector_test_4()
    {
    aligned_vector<double> v(5, 3.14);
    TEST_EQ(5, v.size());
    TEST_EQ(3.14, v[0]);
    TEST_EQ(3.14, v[1]);
    TEST_EQ(3.14, v[2]);
    TEST_EQ(3.14, v[3]);
    TEST_EQ(3.14, v[4]);
    }

  void run_all_vector_tests()
    {
    aligned_vector_test_1();
    aligned_vector_test_2();
    aligned_vector_test_3();
    aligned_vector_test_4();
    }

  void sphere_intersection_1()
    {
    vec3<float4> sphere_origin(float4(3.0), float4(0.0), float4(0.0));
    float4 sphere_radius(1.0);
    vec3<float4> ray_orig(float4(0.0), float4(0.0), float4(0.0));
    vec3<float4> ray_dir(float4(1.0), float4(0.0), float4(0.0));
    float4 t_near(0.0);
    float4 t_far(5.0);
    spherehit4 h;
    intersect_sphere(sphere_origin, sphere_radius, ray_orig, ray_dir, t_near, t_far, h);
    TEST_ASSERT(all(h.found));
    for (int i = 0; i < 4; ++i)
      TEST_EQ_CLOSE(h.distance[i], 2.f, 1e-5);
    }

  void sphere_intersection_2()
    {
    vec3<float4> sphere_origin(float4(-0.5), float4(0.0), float4(0.0));
    float4 sphere_radius(1.0);
    vec3<float4> ray_orig(float4(0.0), float4(0.0), float4(0.0));
    vec3<float4> ray_dir(float4(1.0), float4(0.0), float4(0.0));
    float4 t_near(0.0);
    float4 t_far(5.0);
    spherehit4 h;
    intersect_sphere(sphere_origin, sphere_radius, ray_orig, ray_dir, t_near, t_far, h);
    TEST_ASSERT(all(h.found));
    for (int i = 0; i < 4; ++i)
      TEST_EQ(h.distance[i], 0.5f);
    }

  void sphere_intersection_3()
    {
    vec3<float4> sphere_origin(float4(3.0), float4(0.0), float4(0.0));
    float4 sphere_radius(1.0, 2.0, 2.5, 0.5);
    vec3<float4> ray_orig(float4(0.0), float4(0.0), float4(0.0));
    vec3<float4> ray_dir(float4(1.0), float4(0.0), float4(0.0));
    float4 t_near(0.0);
    float4 t_far(5.0);
    spherehit4 h;
    intersect_sphere(sphere_origin, sphere_radius, ray_orig, ray_dir, t_near, t_far, h);
    TEST_ASSERT(all(h.found));
    TEST_EQ(h.distance[0], 2.f);
    TEST_EQ(h.distance[1], 1.f);
    TEST_EQ(h.distance[2], 0.5f);
    TEST_EQ(h.distance[3], 2.5f);
    }

  void sphere_intersection_4()
    {
    vec3<float4> sphere_origin(float4(-1.5f), float4(0.0), float4(0.0));
    float4 sphere_radius(1.0, 2.0, 2.5, 0.5);
    vec3<float4> ray_orig(float4(0.0), float4(0.0), float4(0.0));
    vec3<float4> ray_dir(float4(1.0), float4(0.0), float4(0.0));
    float4 t_near(0.0);
    float4 t_far(5.0);
    spherehit4 h;
    intersect_sphere(sphere_origin, sphere_radius, ray_orig, ray_dir, t_near, t_far, h);
    TEST_EQ(0, h.found[0]);
    TEST_EQ(-1, h.found[1]);
    TEST_EQ(-1, h.found[2]);
    TEST_EQ(0, h.found[3]);
    TEST_EQ(h.distance[1], 0.5f);
    TEST_EQ(h.distance[2], 1.f);    
    }

  void run_all_sphere_intersection_tests()
    {
    sphere_intersection_1();
    sphere_intersection_2();
    sphere_intersection_3();
    sphere_intersection_4();
    }

  void run_all_quaternion_tests()
    {
    std::vector<float> rx, ry, rz;
    rx = {{0.0, 0.1, 0.2, 0.3, 0.4, -0.1, -0.2, -0.3, -0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -0.5, -0.6, -0.7, -0.8, -0.9, 1.0}};
    ry = { {0.0, 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, -0.1, -0.1, -0.2, -0.2, -0.3, -0.3, -0.4, -0.4, 0.7, 0.8, 0.9, 1.0} };
    rz = { {1.0, 0.9, 0.8, 0.7, -1.0, -0.9, -0.8, -0.7, 0.6, 0.5, 0.4, 0.3, -0.6, -0.5, -0.4, -0.3, 0.1, 0.0, -0.1, -0.9} };
    for (size_t i = 0; i < rx.size(); ++i)
      {
      float x = rx[i];
      float y = ry[i];
      float z = rz[i];
      jtk::float4x4 rot = jtk::compute_from_roll_pitch_yaw_transformation(x, y, z, 0.f, 0.f, 0.f);
      float tx,ty,tz;
      jtk::solve_roll_pitch_yaw_transformation(x, y, z, tx, ty, tz, rot);
      TEST_EQ_CLOSE(x, rx[i], 1e-5f);
      TEST_EQ_CLOSE(y, ry[i], 1e-5f);
      TEST_EQ_CLOSE(z, rz[i], 1e-5f);
      jtk::float4 q = jtk::rotation_to_quaternion(rot);
      jtk::float4 q2 = jtk::roll_pitch_yaw_to_quaternion(x, y, z);
      for (int j = 0; j < 4; ++j)
        TEST_EQ_CLOSE(q[j], q2[j], 1e-5f);
      jtk::float4x4 rot2 = jtk::quaternion_to_rotation(q);
      for (int j = 0; j < 16; ++j)
        TEST_EQ_CLOSE(rot[j], rot2[j], 1e-5f);
      }
    }
  }


void run_all_qbvh_tests()
  {
  using namespace jtk;
  run_all_quaternion_tests();
  run_all_simd_tests();
  run_all_bvh_tests();
  run_all_parallel_partition_tests();
  run_all_vector_tests();
  run_all_sphere_intersection_tests();
  }
