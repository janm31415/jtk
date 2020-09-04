///////////////////////////////////////////////////////////////////////////////
// Includes
/////////////////////////////////////////////////////////////////////////////////

#include "mat_tests.h"

#include <stdint.h>
#include <sstream>
#include <iomanip>

#include "../jtk/mat.h"

#include "test_assert.h"

namespace jtk
  {

  namespace
    {
    struct matrix_fixture
      {
      matrix_fixture() : m(2, 3), m16(2, 3)
        {
        float v = 0.f;
        for (auto& value : m)
          {
          value = v;
          v += 1.f;
          }
        v = 0.f;
        for (auto& value : m16)
          {
          value = v;
          v += 1.f;
          }
        }

      matf m;
      matf16 m16;

      };

    mat3 make_vec3(double a, double b, double c)
      {
      mat3 v(3);
      v(0) = a;
      v(1) = b;
      v(2) = c;
      return v;
      }

    struct matrix_construction_vector : public matrix_fixture {
      void test()
        {
        TEST_EQ(2, m.rows());
        TEST_EQ(3, m.cols());
        TEST_EQ(0.f, m(0, 0));
        TEST_EQ(1.f, m(0, 1));
        TEST_EQ(2.f, m(0, 2));
        TEST_EQ(3.f, m(1, 0));
        TEST_EQ(4.f, m(1, 1));
        TEST_EQ(5.f, m(1, 2));
        TEST_EQ(0.f, m[0][0]);
        TEST_EQ(1.f, m[0][1]);
        TEST_EQ(2.f, m[0][2]);
        TEST_EQ(3.f, m[1][0]);
        TEST_EQ(4.f, m[1][1]);
        TEST_EQ(5.f, m[1][2]);
        TEST_ASSERT(!m.empty());
        }
      };

    struct matrix_construction_array : public matrix_fixture {
      void test()
        {
        TEST_EQ(2, m16.rows());
        TEST_EQ(3, m16.cols());
        TEST_EQ(0.f, m16(0, 0));
        TEST_EQ(1.f, m16(0, 1));
        TEST_EQ(2.f, m16(0, 2));
        TEST_EQ(3.f, m16(1, 0));
        TEST_EQ(4.f, m16(1, 1));
        TEST_EQ(5.f, m16(1, 2));
        TEST_EQ(0.f, m16[0][0]);
        TEST_EQ(1.f, m16[0][1]);
        TEST_EQ(2.f, m16[0][2]);
        TEST_EQ(3.f, m16[1][0]);
        TEST_EQ(4.f, m16[1][1]);
        TEST_EQ(5.f, m16[1][2]);
        TEST_ASSERT(!m16.empty());
        }
      };

    struct copy_constructor : public matrix_fixture {
      void test()
        {
        matf nm(m);
        TEST_EQ(2, nm.rows());
        TEST_EQ(3, nm.cols());
        TEST_EQ(0.f, nm(0, 0));
        TEST_EQ(1.f, nm(0, 1));
        TEST_EQ(2.f, nm(0, 2));
        TEST_EQ(3.f, nm(1, 0));
        TEST_EQ(4.f, nm(1, 1));
        TEST_EQ(5.f, nm(1, 2));
        TEST_EQ(0.f, nm[0][0]);
        TEST_EQ(1.f, nm[0][1]);
        TEST_EQ(2.f, nm[0][2]);
        TEST_EQ(3.f, nm[1][0]);
        TEST_EQ(4.f, nm[1][1]);
        TEST_EQ(5.f, nm[1][2]);
        TEST_ASSERT(!nm.empty());
        }
      };

    struct copy_constructor_array : public matrix_fixture {
      void test()
        {
        matf16 nm(m16);
        TEST_EQ(2, nm.rows());
        TEST_EQ(3, nm.cols());
        TEST_EQ(0.f, nm(0, 0));
        TEST_EQ(1.f, nm(0, 1));
        TEST_EQ(2.f, nm(0, 2));
        TEST_EQ(3.f, nm(1, 0));
        TEST_EQ(4.f, nm(1, 1));
        TEST_EQ(5.f, nm(1, 2));
        TEST_EQ(0.f, nm[0][0]);
        TEST_EQ(1.f, nm[0][1]);
        TEST_EQ(2.f, nm[0][2]);
        TEST_EQ(3.f, nm[1][0]);
        TEST_EQ(4.f, nm[1][1]);
        TEST_EQ(5.f, nm[1][2]);
        TEST_ASSERT(!nm.empty());
        }
      };

    struct copy_assignment_constructor : public matrix_fixture {
      void test()
        {
        matf nm;
        TEST_ASSERT(nm.empty());
        nm = m;
        TEST_EQ(2, nm.rows());
        TEST_EQ(3, nm.cols());
        TEST_EQ(0.f, nm(0, 0));
        TEST_EQ(1.f, nm(0, 1));
        TEST_EQ(2.f, nm(0, 2));
        TEST_EQ(3.f, nm(1, 0));
        TEST_EQ(4.f, nm(1, 1));
        TEST_EQ(5.f, nm(1, 2));
        TEST_EQ(0.f, nm[0][0]);
        TEST_EQ(1.f, nm[0][1]);
        TEST_EQ(2.f, nm[0][2]);
        TEST_EQ(3.f, nm[1][0]);
        TEST_EQ(4.f, nm[1][1]);
        TEST_EQ(5.f, nm[1][2]);
        TEST_ASSERT(!nm.empty());
        }
      };

    struct copy_assignment_array : public matrix_fixture {
      void test()
        {
        matf16 nm;
        TEST_ASSERT(nm.empty());
        nm = m16;
        TEST_EQ(2, nm.rows());
        TEST_EQ(3, nm.cols());
        TEST_EQ(0.f, nm(0, 0));
        TEST_EQ(1.f, nm(0, 1));
        TEST_EQ(2.f, nm(0, 2));
        TEST_EQ(3.f, nm(1, 0));
        TEST_EQ(4.f, nm(1, 1));
        TEST_EQ(5.f, nm(1, 2));
        TEST_EQ(0.f, nm[0][0]);
        TEST_EQ(1.f, nm[0][1]);
        TEST_EQ(2.f, nm[0][2]);
        TEST_EQ(3.f, nm[1][0]);
        TEST_EQ(4.f, nm[1][1]);
        TEST_EQ(5.f, nm[1][2]);
        TEST_ASSERT(!nm.empty());
        }
      };

    void vector_test()
      {
      mat v(10);
      double d = 0.0;
      for (auto& value : v)
        {
        value = d;
        d += 1.0;
        }
      for (int i = 0; i < 10; ++i)
        TEST_EQ((double)i, v(i));
      }

    void vector_test_array()
      {
      mat10 v(10);
      double d = 0.0;
      for (auto& value : v)
        {
        value = d;
        d += 1.0;
        }
      for (int i = 0; i < 10; ++i)
        TEST_EQ((double)i, v(i));
      }


    void svd_1()
      {
      auto p0 = make_vec3(1.0, 0.0, 0.0);
      auto p1 = make_vec3(0.0, 1.0, 0.0);
      auto p2 = make_vec3(0.0, 0.0, 1.0);
      std::vector<mat3> pts;
      pts.push_back(p0);
      pts.push_back(p1);
      pts.push_back(p2);
      auto c = make_vec3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
      mat a(3, 3);
      for (int i = 0; i < 3; ++i)
        {
        for (int j = 0; j < 3; ++j)
          a[i][j] = 0.0;
        }
      for (size_t k = 0; k < pts.size(); ++k)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            a[i][j] += (pts[k](i) - c(i))*(pts[k](j) - c(j));
      mat w, v;
      bool res = svd(a, w, v);
      TEST_ASSERT(res);

      TEST_EQ_CLOSE(1.0, w(0), 1e-12);
      TEST_EQ_CLOSE(1.0, w(1), 1e-12);
      TEST_EQ_CLOSE(0.0, w(2), 1e-12);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ_CLOSE(a[i][j], v[i][j], 1e-12);
      }


    void svd_2()
      {
      int m = 4;
      int n = 3;
      matf A(m, n);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = float(i);
      matf A_copy(A);

      matf sigma;
      matf V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          float val = 0.f;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-5);
          }
        }
      }

    void svd_3()
      {
      int m = 3;
      int n = 4;
      matf A(m, n);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) % 3][(i - 1) / 3] = float(i);
      matf A_copy(A);

      matf sigma;
      matf V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          float val = 0.f;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-5);
          }
        }
      }

    void pseudo_inverse_1()
      {
      int m = 3;
      int n = 3;
      mat A(m, n);
      for (int i = 1; i <= 9; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = std::sqrt(double(i - 1));
      mat A_copy(A);

      mat inv;

      int p = pseudo_inverse(inv, A, 1e-5);
      TEST_EQ(3, p);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < m; ++c)
          {
          double val = 0.0;
          for (int k = 0; k < n; ++k)
            {
            val += A_copy[r][k] * inv[k][c];
            }
          TEST_EQ_CLOSE(r == c ? 1.0 : 0.0, val, 1e-12);
          }
        }
      }


    void pseudo_inverse_2()
      {
      int m = 8;
      int n = 6;
      mat A(m, n);

      A[0][0] = 64;
      A[0][1] = 2;
      A[0][2] = 3;
      A[0][3] = 61;
      A[0][4] = 60;
      A[0][5] = 6;
      A[1][0] = 9;
      A[1][1] = 55;
      A[1][2] = 54;
      A[1][3] = 12;
      A[1][4] = 13;
      A[1][5] = 51;
      A[2][0] = 17;
      A[2][1] = 47;
      A[2][2] = 46;
      A[2][3] = 20;
      A[2][4] = 21;
      A[2][5] = 43;
      A[3][0] = 40;
      A[3][1] = 26;
      A[3][2] = 27;
      A[3][3] = 37;
      A[3][4] = 36;
      A[3][5] = 30;
      A[4][0] = 32;
      A[4][1] = 34;
      A[4][2] = 35;
      A[4][3] = 29;
      A[4][4] = 28;
      A[4][5] = 38;
      A[5][0] = 41;
      A[5][1] = 23;
      A[5][2] = 22;
      A[5][3] = 44;
      A[5][4] = 45;
      A[5][5] = 19;
      A[6][0] = 49;
      A[6][1] = 15;
      A[6][2] = 14;
      A[6][3] = 52;
      A[6][4] = 53;
      A[6][5] = 11;
      A[7][0] = 8;
      A[7][1] = 58;
      A[7][2] = 59;
      A[7][3] = 5;
      A[7][4] = 4;
      A[7][5] = 62;


      mat A_copy(A);

      mat sigma;
      mat V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          double val = 0.0;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-12);
          }
        }

      A = A_copy;
      mat inv;
      int p = pseudo_inverse(inv, A, 1e-10);
      TEST_EQ(3, p);

      double b[8] = { 260.0, 260.0, 260.0, 260.0, 260.0, 260.0, 260.0, 260.0 };

      double x[6];
      for (int j = 0; j < 6; ++j)
        {
        x[j] = 0.0;
        for (int i = 0; i < 8; ++i)
          x[j] += inv[j][i] * b[i];
        }

      double n2 = 0.0;
      double error[8];
      for (int i = 0; i < 8; ++i)
        {
        error[i] = -b[i];
        for (int j = 0; j < 6; ++j)
          error[i] += A_copy[i][j] * x[j];
        n2 += error[i] * error[i];
        }
      TEST_ASSERT(n2 < 1e-8);
      }

    void pseudo_inverse_3()
      {
      int m = 6;
      int n = 8;
      mat A(m, n);

      A[0][0] = 64;
      A[1][0] = 2;
      A[2][0] = 3;
      A[3][0] = 61;
      A[4][0] = 60;
      A[5][0] = 6;
      A[0][1] = 9;
      A[1][1] = 55;
      A[2][1] = 54;
      A[3][1] = 12;
      A[4][1] = 13;
      A[5][1] = 51;
      A[0][2] = 17;
      A[1][2] = 47;
      A[2][2] = 46;
      A[3][2] = 20;
      A[4][2] = 21;
      A[5][2] = 43;
      A[0][3] = 40;
      A[1][3] = 26;
      A[2][3] = 27;
      A[3][3] = 37;
      A[4][3] = 36;
      A[5][3] = 30;
      A[0][4] = 32;
      A[1][4] = 34;
      A[2][4] = 35;
      A[3][4] = 29;
      A[4][4] = 28;
      A[5][4] = 38;
      A[0][5] = 41;
      A[1][5] = 23;
      A[2][5] = 22;
      A[3][5] = 44;
      A[4][5] = 45;
      A[5][5] = 19;
      A[0][6] = 49;
      A[1][6] = 15;
      A[2][6] = 14;
      A[3][6] = 52;
      A[4][6] = 53;
      A[5][6] = 11;
      A[0][7] = 8;
      A[1][7] = 58;
      A[2][7] = 59;
      A[3][7] = 5;
      A[4][7] = 4;
      A[5][7] = 62;


      mat A_copy(A);

      mat sigma;
      mat V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          double val = 0.0;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-12);
          }
        }

      A = A_copy;
      mat inv;
      int p = pseudo_inverse(inv, A, 1e-10);
      TEST_EQ(3, p);
      double b[6] = { 260.0, 260.0, 260.0, 260.0, 260.0, 260.0 };

      double x[8];
      for (int j = 0; j < 8; ++j)
        {
        x[j] = 0.0;
        for (int i = 0; i < 6; ++i)
          x[j] += inv[j][i] * b[i];
        }

      double n2 = 0.0;
      double error[6];
      for (int i = 0; i < 6; ++i)
        {
        error[i] = -b[i];
        for (int j = 0; j < 8; ++j)
          error[i] += A_copy[i][j] * x[j];
        n2 += error[i] * error[i];
        }
      TEST_ASSERT(n2 < 1e-8);
      }

    void svd_1_array()
      {
      auto p0 = make_vec3(1.0, 0.0, 0.0);
      auto p1 = make_vec3(0.0, 1.0, 0.0);
      auto p2 = make_vec3(0.0, 0.0, 1.0);
      std::vector<mat3> pts;
      pts.push_back(p0);
      pts.push_back(p1);
      pts.push_back(p2);
      auto c = make_vec3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
      mat9 a(3, 3);
      for (int i = 0; i < 3; ++i)
        {
        for (int j = 0; j < 3; ++j)
          a[i][j] = 0.0;
        }
      for (size_t k = 0; k < pts.size(); ++k)
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            a[i][j] += (pts[k](i) - c(i))*(pts[k](j) - c(j));
      mat3 w;
      mat9 v;
      bool res = svd(a, w, v);
      TEST_ASSERT(res);

      TEST_EQ_CLOSE(1.0, w(0), 1e-8);
      TEST_EQ_CLOSE(1.0, w(1), 1e-8);
      TEST_EQ_CLOSE(0.0, w(2), 1e-8);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ_CLOSE(a[i][j], v[i][j], 1e-8);
      }


    void svd_2_array()
      {
      int m = 4;
      int n = 3;
      matf12 A(m, n);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = float(i);
      matf12 A_copy(A);

      matf3 sigma;
      matf9 V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          float val = 0.f;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-5);
          }
        }
      }

    void svd_3_array()
      {
      int m = 3;
      int n = 4;
      matf12 A(m, n);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) % 3][(i - 1) / 3] = float(i);
      matf12 A_copy(A);

      matf4 sigma;
      matf16 V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          float val = 0.f;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-5);
          }
        }
      }

    void pseudo_inverse_1_array()
      {
      int m = 3;
      int n = 3;
      mat9 A(m, n);
      for (int i = 1; i <= 9; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = std::sqrt(double(i - 1));
      mat9 A_copy(A);

      mat9 inv;

      int res = pseudo_inverse(inv, A, 1e-5);
      TEST_EQ(3, res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < m; ++c)
          {
          double val = 0.0;
          for (int k = 0; k < n; ++k)
            {
            val += A_copy[r][k] * inv[k][c];
            }
          TEST_EQ_CLOSE(r == c ? 1.0 : 0.0, val, 1e-12);
          }
        }
      }

    void pseudo_inverse_2_array()
      {
      int m = 8;
      int n = 6;
      mat48 A(m, n);

      A[0][0] = 64;
      A[0][1] = 2;
      A[0][2] = 3;
      A[0][3] = 61;
      A[0][4] = 60;
      A[0][5] = 6;
      A[1][0] = 9;
      A[1][1] = 55;
      A[1][2] = 54;
      A[1][3] = 12;
      A[1][4] = 13;
      A[1][5] = 51;
      A[2][0] = 17;
      A[2][1] = 47;
      A[2][2] = 46;
      A[2][3] = 20;
      A[2][4] = 21;
      A[2][5] = 43;
      A[3][0] = 40;
      A[3][1] = 26;
      A[3][2] = 27;
      A[3][3] = 37;
      A[3][4] = 36;
      A[3][5] = 30;
      A[4][0] = 32;
      A[4][1] = 34;
      A[4][2] = 35;
      A[4][3] = 29;
      A[4][4] = 28;
      A[4][5] = 38;
      A[5][0] = 41;
      A[5][1] = 23;
      A[5][2] = 22;
      A[5][3] = 44;
      A[5][4] = 45;
      A[5][5] = 19;
      A[6][0] = 49;
      A[6][1] = 15;
      A[6][2] = 14;
      A[6][3] = 52;
      A[6][4] = 53;
      A[6][5] = 11;
      A[7][0] = 8;
      A[7][1] = 58;
      A[7][2] = 59;
      A[7][3] = 5;
      A[7][4] = 4;
      A[7][5] = 62;


      mat48 A_copy(A);

      mat48 sigma;
      mat48 V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          double val = 0.0;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-12);
          }
        }

      A = A_copy;
      mat48 inv;
      int p = pseudo_inverse(inv, A, 1e-10);
      TEST_EQ(3, p);

      double b[8] = { 260.0, 260.0, 260.0, 260.0, 260.0, 260.0, 260.0, 260.0 };

      double x[6];
      for (int j = 0; j < 6; ++j)
        {
        x[j] = 0.0;
        for (int i = 0; i < 8; ++i)
          x[j] += inv[j][i] * b[i];
        }

      double n2 = 0.0;
      double error[8];
      for (int i = 0; i < 8; ++i)
        {
        error[i] = -b[i];
        for (int j = 0; j < 6; ++j)
          error[i] += A_copy[i][j] * x[j];
        n2 += error[i] * error[i];
        }
      TEST_ASSERT(n2 < 1e-8);
      }

    void pseudo_inverse_3_array()
      {
      int m = 6;
      int n = 8;
      mat64 A(m, n);

      A[0][0] = 64;
      A[1][0] = 2;
      A[2][0] = 3;
      A[3][0] = 61;
      A[4][0] = 60;
      A[5][0] = 6;
      A[0][1] = 9;
      A[1][1] = 55;
      A[2][1] = 54;
      A[3][1] = 12;
      A[4][1] = 13;
      A[5][1] = 51;
      A[0][2] = 17;
      A[1][2] = 47;
      A[2][2] = 46;
      A[3][2] = 20;
      A[4][2] = 21;
      A[5][2] = 43;
      A[0][3] = 40;
      A[1][3] = 26;
      A[2][3] = 27;
      A[3][3] = 37;
      A[4][3] = 36;
      A[5][3] = 30;
      A[0][4] = 32;
      A[1][4] = 34;
      A[2][4] = 35;
      A[3][4] = 29;
      A[4][4] = 28;
      A[5][4] = 38;
      A[0][5] = 41;
      A[1][5] = 23;
      A[2][5] = 22;
      A[3][5] = 44;
      A[4][5] = 45;
      A[5][5] = 19;
      A[0][6] = 49;
      A[1][6] = 15;
      A[2][6] = 14;
      A[3][6] = 52;
      A[4][6] = 53;
      A[5][6] = 11;
      A[0][7] = 8;
      A[1][7] = 58;
      A[2][7] = 59;
      A[3][7] = 5;
      A[4][7] = 4;
      A[5][7] = 62;


      mat64 A_copy(A);

      mat64 sigma;
      mat64 V;

      bool res = svd(A, sigma, V);
      TEST_ASSERT(res);

      for (int r = 0; r < m; ++r)
        {
        for (int c = 0; c < n; ++c)
          {
          double val = 0.0;
          for (int k = 0; k < n; ++k)
            {
            val += A[r][k] * sigma(k) * V[c][k];
            }
          TEST_EQ_CLOSE(A_copy[r][c], val, 1e-12);
          }
        }

      A = A_copy;
      mat64 inv;
      int p = pseudo_inverse(inv, A, 1e-10);
      TEST_EQ(3, p);
      double b[6] = { 260.0, 260.0, 260.0, 260.0, 260.0, 260.0 };

      double x[8];
      for (int j = 0; j < 8; ++j)
        {
        x[j] = 0.0;
        for (int i = 0; i < 6; ++i)
          x[j] += inv[j][i] * b[i];
        }

      double n2 = 0.0;
      double error[6];
      for (int i = 0; i < 6; ++i)
        {
        error[i] = -b[i];
        for (int j = 0; j < 8; ++j)
          error[i] += A_copy[i][j] * x[j];
        n2 += error[i] * error[i];
        }
      TEST_ASSERT(n2 < 1e-8);
      }

    struct adding : public matrix_fixture {
      void test()
        {
        matf res = m + m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        }
      };

    struct adding2 : public matrix_fixture {
      void test()
        {
        matf res(m + m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        }
      };

    struct adding3 : public matrix_fixture {
      void test()
        {
        auto expr = m + m;
        matf res;
        res.assign(expr);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        }
      };

    struct adding4 : public matrix_fixture {
      void test()
        {
        auto expr = m + m;
        matf res = m + expr;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(6.f, res(0, 2));
        TEST_EQ(9.f, res(1, 0));
        TEST_EQ(12.f, res(1, 1));
        TEST_EQ(15.f, res(1, 2));
        }
      };

    struct adding5 : public matrix_fixture {
      void test()
        {
        auto expr = m + m;
        matf res = expr + m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(6.f, res(0, 2));
        TEST_EQ(9.f, res(1, 0));
        TEST_EQ(12.f, res(1, 1));
        TEST_EQ(15.f, res(1, 2));
        }
      };

    struct adding6 : public matrix_fixture {
      void test()
        {
        auto expr1 = m + m;
        auto expr2 = m + m;
        matf res = expr1 + expr2;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        res = expr1 + expr1;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        }
      };


    struct subtracting1 : public matrix_fixture {
      void test()
        {
        auto expr = m - m;
        matf res;
        res.assign(expr);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(0.f, res(0, 1));
        TEST_EQ(0.f, res(0, 2));
        TEST_EQ(0.f, res(1, 0));
        TEST_EQ(0.f, res(1, 1));
        TEST_EQ(0.f, res(1, 2));
        }
      };

    struct subtracting2 : public matrix_fixture {
      void test()
        {
        auto expr = m - m;
        matf res = m - expr;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(1.f, res(0, 1));
        TEST_EQ(2.f, res(0, 2));
        TEST_EQ(3.f, res(1, 0));
        TEST_EQ(4.f, res(1, 1));
        TEST_EQ(5.f, res(1, 2));
        }
      };

    struct subtracting3 : public matrix_fixture {
      void test()
        {
        auto expr = m - m;
        matf res = expr - m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(0, 1));
        TEST_EQ(-2.f, res(0, 2));
        TEST_EQ(-3.f, res(1, 0));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(1, 2));
        }
      };

    struct subtracting4 : public matrix_fixture {
      void test()
        {
        auto expr1 = m - m;
        auto expr2 = m - m;
        matf res = expr1 - expr2;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(0.f, res(0, 1));
        TEST_EQ(0.f, res(0, 2));
        TEST_EQ(0.f, res(1, 0));
        TEST_EQ(0.f, res(1, 1));
        TEST_EQ(0.f, res(1, 2));
        }
      };

    struct negating1 : public matrix_fixture {
      void test()
        {
        matf res = -m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(0, 1));
        TEST_EQ(-2.f, res(0, 2));
        TEST_EQ(-3.f, res(1, 0));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(1, 2));
        }
      };

    struct negating2 : public matrix_fixture {
      void test()
        {
        auto expr = m + m;
        matf res = -expr;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-2.f, res(0, 1));
        TEST_EQ(-4.f, res(0, 2));
        TEST_EQ(-6.f, res(1, 0));
        TEST_EQ(-8.f, res(1, 1));
        TEST_EQ(-10.f, res(1, 2));
        }
      };

    struct negating3 : public matrix_fixture {
      void test()
        {
        matf res = -transpose(m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(1, 0));
        TEST_EQ(-2.f, res(2, 0));
        TEST_EQ(-3.f, res(0, 1));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(2, 1));
        }
      };

    void transpose_aliasing_problem()
      {
      jtk::mat m(3, 3);
      m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
      m = jtk::transpose(m);
      TEST_EQ(1.0, m(0, 0));
      TEST_EQ(2.0, m(1, 0));
      TEST_EQ(3.0, m(2, 0));
      TEST_EQ(4.0, m(0, 1));
      TEST_EQ(5.0, m(1, 1));
      TEST_EQ(6.0, m(2, 1));
      TEST_EQ(7.0, m(0, 2));
      TEST_EQ(8.0, m(1, 2));
      TEST_EQ(9.0, m(2, 2));
      }

    void transpose_aliasing_problem_2()
      {
      jtk::mat m(3, 3);
      m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
      m = -jtk::transpose(m);
      TEST_EQ(-1.0, m(0, 0));
      TEST_EQ(-2.0, m(1, 0));
      TEST_EQ(-3.0, m(2, 0));
      TEST_EQ(-4.0, m(0, 1));
      TEST_EQ(-5.0, m(1, 1));
      TEST_EQ(-6.0, m(2, 1));
      TEST_EQ(-7.0, m(0, 2));
      TEST_EQ(-8.0, m(1, 2));
      TEST_EQ(-9.0, m(2, 2));
      }

    struct scalarmultiply : public matrix_fixture {
      void test()
        {
        matf res = m * 2.f;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        res = 3.f * m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(6.f, res(0, 2));
        TEST_EQ(9.f, res(1, 0));
        TEST_EQ(12.f, res(1, 1));
        TEST_EQ(15.f, res(1, 2));
        res = (m + m) * 2.f;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        res = 3.f * (m + m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(6.f, res(0, 1));
        TEST_EQ(12.f, res(0, 2));
        TEST_EQ(18.f, res(1, 0));
        TEST_EQ(24.f, res(1, 1));
        TEST_EQ(30.f, res(1, 2));
        }
      };

    void scalardivision()
      {
      matf m(3, 1);
      m(0) = 8.f;
      m(1) = 4.f;
      m(2) = 2.f;
      matf res = m / 2.f;
      TEST_EQ(4.f, res(0));
      TEST_EQ(2.f, res(1));
      TEST_EQ(1.f, res(2));
      res = 8.f / m;
      TEST_EQ(1.f, res(0));
      TEST_EQ(2.f, res(1));
      TEST_EQ(4.f, res(2));
      res = (m + m) / 2.f;
      TEST_EQ(8.f, res(0));
      TEST_EQ(4.f, res(1));
      TEST_EQ(2.f, res(2));
      res = 8.f / (m + m);
      TEST_EQ(0.5f, res(0));
      TEST_EQ(1.f, res(1));
      TEST_EQ(2.f, res(2));
      }

    void scalaraddition()
      {
      matf m(3, 1);
      m(0) = 8.f;
      m(1) = 4.f;
      m(2) = 2.f;
      matf res = m + 1.f;
      TEST_EQ(9.f, res(0));
      TEST_EQ(5.f, res(1));
      TEST_EQ(3.f, res(2));
      res = 1.f + m;
      TEST_EQ(9.f, res(0));
      TEST_EQ(5.f, res(1));
      TEST_EQ(3.f, res(2));
      res = (m + m) + 1.f;
      TEST_EQ(17.f, res(0));
      TEST_EQ(9.f, res(1));
      TEST_EQ(5.f, res(2));
      res = 1.f + (m + m);
      TEST_EQ(17.f, res(0));
      TEST_EQ(9.f, res(1));
      TEST_EQ(5.f, res(2));
      }

    void scalarsubtraction()
      {
      matf m(3, 1);
      m(0) = 8.f;
      m(1) = 4.f;
      m(2) = 2.f;
      matf res = m - 1.f;
      TEST_EQ(7.f, res(0));
      TEST_EQ(3.f, res(1));
      TEST_EQ(1.f, res(2));
      res = 1.f - m;
      TEST_EQ(-7.f, res(0));
      TEST_EQ(-3.f, res(1));
      TEST_EQ(-1.f, res(2));
      res = (m + m) - 1.f;
      TEST_EQ(15.f, res(0));
      TEST_EQ(7.f, res(1));
      TEST_EQ(3.f, res(2));
      res = 1.f - (m + m);
      TEST_EQ(-15.f, res(0));
      TEST_EQ(-7.f, res(1));
      TEST_EQ(-3.f, res(2));
      }

    void matrix_multiplication()
      {
      mat m1(2, 3);
      mat m2(3, 2);
      mat m1zero(2, 3);
      mat m2zero(3, 2);

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      mat m = m1 * m2;
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = m1 * (m2 + m2zero);
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = m1 * (m2 + m2zero + m2zero);
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = (m1 + m1zero) * m2;
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = (m1 + m1zero + m1zero) * m2;
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = (m1 + m1zero + m1zero) * (m2 + 1.0 - 1.0);
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));
      }

    struct adding_array : public matrix_fixture {
      void test()
        {
        matf25 res = m16 + m16;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        }
      };

    struct adding2_array : public matrix_fixture {
      void test()
        {
        matf25 res(m16 + m16);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        }
      };

    struct adding3_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 + m16;
        matf25 res;
        res.assign(expr);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        }
      };

    struct adding4_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 + m16;
        matf25 res = m16 + expr;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(6.f, res(0, 2));
        TEST_EQ(9.f, res(1, 0));
        TEST_EQ(12.f, res(1, 1));
        TEST_EQ(15.f, res(1, 2));
        }
      };

    struct adding5_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 + m16;
        matf25 res = expr + m16;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(6.f, res(0, 2));
        TEST_EQ(9.f, res(1, 0));
        TEST_EQ(12.f, res(1, 1));
        TEST_EQ(15.f, res(1, 2));
        }
      };

    struct adding6_array : public matrix_fixture {
      void test()
        {
        auto expr1 = m16 + m16;
        auto expr2 = m16 + m16;
        matf25 res = expr1 + expr2;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        res = expr1 + expr1;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        }
      };

    struct subtracting1_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 - m16;
        matf25 res;
        res.assign(expr);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(0.f, res(0, 1));
        TEST_EQ(0.f, res(0, 2));
        TEST_EQ(0.f, res(1, 0));
        TEST_EQ(0.f, res(1, 1));
        TEST_EQ(0.f, res(1, 2));
        }
      };

    struct subtracting2_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 - m16;
        matf25 res = m16 - expr;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(1.f, res(0, 1));
        TEST_EQ(2.f, res(0, 2));
        TEST_EQ(3.f, res(1, 0));
        TEST_EQ(4.f, res(1, 1));
        TEST_EQ(5.f, res(1, 2));
        }
      };

    struct subtracting3_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 - m16;
        matf25 res = expr - m16;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(0, 1));
        TEST_EQ(-2.f, res(0, 2));
        TEST_EQ(-3.f, res(1, 0));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(1, 2));
        }
      };
    struct subtracting4_array : public matrix_fixture {
      void test()
        {
        auto expr1 = m16 - m16;
        auto expr2 = m16 - m16;
        matf25 res = expr1 - expr2;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(0.f, res(0, 1));
        TEST_EQ(0.f, res(0, 2));
        TEST_EQ(0.f, res(1, 0));
        TEST_EQ(0.f, res(1, 1));
        TEST_EQ(0.f, res(1, 2));
        }
      };

    struct negating1_array : public matrix_fixture {
      void test()
        {
        matf25 res = -m16;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(0, 1));
        TEST_EQ(-2.f, res(0, 2));
        TEST_EQ(-3.f, res(1, 0));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(1, 2));
        }
      };

    struct negating2_array : public matrix_fixture {
      void test()
        {
        auto expr = m16 + m16;
        matf25 res = -expr;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-2.f, res(0, 1));
        TEST_EQ(-4.f, res(0, 2));
        TEST_EQ(-6.f, res(1, 0));
        TEST_EQ(-8.f, res(1, 1));
        TEST_EQ(-10.f, res(1, 2));
        }
      };

    struct scalarmultiply_array : public matrix_fixture {
      void test()
        {
        matf25 res = m16 * 2.f;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));
        res = 3.f * m16;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(6.f, res(0, 2));
        TEST_EQ(9.f, res(1, 0));
        TEST_EQ(12.f, res(1, 1));
        TEST_EQ(15.f, res(1, 2));
        res = (m16 + m16) * 2.f;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        res = 3.f * (m16 + m16);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(6.f, res(0, 1));
        TEST_EQ(12.f, res(0, 2));
        TEST_EQ(18.f, res(1, 0));
        TEST_EQ(24.f, res(1, 1));
        TEST_EQ(30.f, res(1, 2));
        }
      };
    void scalardivision_array()
      {
      matf25 m(3, 1);
      m(0) = 8.f;
      m(1) = 4.f;
      m(2) = 2.f;
      matf25 res = m / 2.f;
      TEST_EQ(4.f, res(0));
      TEST_EQ(2.f, res(1));
      TEST_EQ(1.f, res(2));
      res = 8.f / m;
      TEST_EQ(1.f, res(0));
      TEST_EQ(2.f, res(1));
      TEST_EQ(4.f, res(2));
      res = (m + m) / 2.f;
      TEST_EQ(8.f, res(0));
      TEST_EQ(4.f, res(1));
      TEST_EQ(2.f, res(2));
      res = 8.f / (m + m);
      TEST_EQ(0.5f, res(0));
      TEST_EQ(1.f, res(1));
      TEST_EQ(2.f, res(2));
      }

    void scalaraddition_array()
      {
      matf25 m(3, 1);
      m(0) = 8.f;
      m(1) = 4.f;
      m(2) = 2.f;
      matf25 res = m + 1.f;
      TEST_EQ(9.f, res(0));
      TEST_EQ(5.f, res(1));
      TEST_EQ(3.f, res(2));
      res = 1.f + m;
      TEST_EQ(9.f, res(0));
      TEST_EQ(5.f, res(1));
      TEST_EQ(3.f, res(2));
      res = (m + m) + 1.f;
      TEST_EQ(17.f, res(0));
      TEST_EQ(9.f, res(1));
      TEST_EQ(5.f, res(2));
      res = 1.f + (m + m);
      TEST_EQ(17.f, res(0));
      TEST_EQ(9.f, res(1));
      TEST_EQ(5.f, res(2));
      }

    void scalarsubtraction_array()
      {
      matf25 m(3, 1);
      m(0) = 8.f;
      m(1) = 4.f;
      m(2) = 2.f;
      matf25 res = m - 1.f;
      TEST_EQ(7.f, res(0));
      TEST_EQ(3.f, res(1));
      TEST_EQ(1.f, res(2));
      res = 1.f - m;
      TEST_EQ(-7.f, res(0));
      TEST_EQ(-3.f, res(1));
      TEST_EQ(-1.f, res(2));
      res = (m + m) - 1.f;
      TEST_EQ(15.f, res(0));
      TEST_EQ(7.f, res(1));
      TEST_EQ(3.f, res(2));
      res = 1.f - (m + m);
      TEST_EQ(-15.f, res(0));
      TEST_EQ(-7.f, res(1));
      TEST_EQ(-3.f, res(2));
      }

    void matrix_multiplication_array()
      {
      mat25 m1(2, 3);
      mat25 m2(3, 2);
      mat25 m1zero = zeros(2, 3);
      mat25 m2zero = zeros(3, 2);

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      mat m = m1 * m2;
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = m1 * (m2 + m2zero);
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = m1 * (m2 + m2zero + m2zero);
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = (m1 + m1zero) * m2;
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = (m1 + m1zero + m1zero) * m2;
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));

      m = (m1 + m1zero + m1zero) * (m2 + 1.0 - 1.0);
      TEST_EQ(2, m.rows());
      TEST_EQ(2, m.cols());
      TEST_EQ(13.0, m(0, 0));
      TEST_EQ(6.0, m(1, 0));
      TEST_EQ(69.0, m(0, 1));
      TEST_EQ(-1.0, m(1, 1));
      }

    struct matrixconversion : public matrix_fixture {
      void test()
        {
        mat9 m9 = m;
        TEST_EQ(2, m9.rows());
        TEST_EQ(3, m9.cols());
        TEST_EQ(0.0, m9(0, 0));
        TEST_EQ(1.0, m9(0, 1));
        TEST_EQ(2.0, m9(0, 2));
        TEST_EQ(3.0, m9(1, 0));
        TEST_EQ(4.0, m9(1, 1));
        TEST_EQ(5.0, m9(1, 2));
        }
      };

    struct matrixtranspose : public matrix_fixture {
      void test()
        {
        matf res = transpose(m);
        TEST_EQ(3, res.rows());
        TEST_EQ(2, res.cols());
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(1.f, res(1, 0));
        TEST_EQ(2.f, res(2, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(4.f, res(1, 1));
        TEST_EQ(5.f, res(2, 1));
        res = transpose(m + m);
        TEST_EQ(3, res.rows());
        TEST_EQ(2, res.cols());
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(1, 0));
        TEST_EQ(4.f, res(2, 0));
        TEST_EQ(6.f, res(0, 1));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(2, 1));
        }
      };

    struct matrixtranspose_array : public matrix_fixture {
      void test()
        {
        matf6 res = transpose(m16);
        TEST_EQ(3, res.rows());
        TEST_EQ(2, res.cols());
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(1.f, res(1, 0));
        TEST_EQ(2.f, res(2, 0));
        TEST_EQ(3.f, res(0, 1));
        TEST_EQ(4.f, res(1, 1));
        TEST_EQ(5.f, res(2, 1));
        res = transpose(m16 + m16);
        TEST_EQ(3, res.rows());
        TEST_EQ(2, res.cols());
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(1, 0));
        TEST_EQ(4.f, res(2, 0));
        TEST_EQ(6.f, res(0, 1));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(2, 1));
        }
      };

    void lsd_test()
      {
      mat A(4, 3);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = double(i);
      mat b(4);
      for (int i = 0; i < 4; ++i)
        b(i) = double(3 * i + 2);
      lsd(A, b, 1e-12);
      TEST_EQ(3, b.rows());
      TEST_EQ(1, b.cols());
      TEST_EQ_CLOSE(1.0 / 3.0, b(0), 1e-12);
      TEST_EQ_CLOSE(1.0 / 3.0, b(1), 1e-12);
      TEST_EQ_CLOSE(1.0 / 3.0, b(2), 1e-12);
      }

    void lsd_test_array()
      {
      mat12 A(4, 3);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = double(i);
      mat4 b(4);
      for (int i = 0; i < 4; ++i)
        b(i) = double(3 * i + 2);
      lsd(A, b, 1e-12);
      TEST_EQ(3, b.rows());
      TEST_EQ(1, b.cols());
      TEST_EQ_CLOSE(1.0 / 3.0, b(0), 1e-12);
      TEST_EQ_CLOSE(1.0 / 3.0, b(1), 1e-12);
      TEST_EQ_CLOSE(1.0 / 3.0, b(2), 1e-12);
      }

    void outputstream()
      {
      mat m(2, 3);
      m(0, 0) = 1;
      m(0, 1) = 10000;
      m(0, 2) = 0.0000143567;
      m(1, 0) = -1.34;
      m(1, 1) = 324.32049823;
      m(1, 2) = 1e-12;
      std::stringstream str;
      str << m;
      TEST_EQ(72, str.str().length());
      str.clear();
      str.str("");
      str << std::scientific;
      str << std::setprecision(20);
      str << m;
      TEST_EQ(168, str.str().length());
      }

    void inputstream()
      {
      mat m(2, 3);
      m << 1, 2, 3, 4, 5, 6;
      TEST_EQ(1.0, m(0, 0));
      TEST_EQ(2.0, m(0, 1));
      TEST_EQ(3.0, m(0, 2));
      TEST_EQ(4.0, m(1, 0));
      TEST_EQ(5.0, m(1, 1));
      TEST_EQ(6.0, m(1, 2));
      }

    struct matrix_add_assignment : public matrix_fixture {
      void test()
        {
        matf res(m);
        res += m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));

        res += (m + m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        }
      };

    struct matrix_subtract_assignment : public matrix_fixture {
      void test()
        {
        matf res(2, 3);
        res -= m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(0, 1));
        TEST_EQ(-2.f, res(0, 2));
        TEST_EQ(-3.f, res(1, 0));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(1, 2));
        res -= (m + m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-3.f, res(0, 1));
        TEST_EQ(-6.f, res(0, 2));
        TEST_EQ(-9.f, res(1, 0));
        TEST_EQ(-12.f, res(1, 1));
        TEST_EQ(-15.f, res(1, 2));
        }
      };

    void matrix_multiply_assignment()
      {
      mat m1(2, 3);
      mat m2(3, 2);

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      m1 *= m2;
      TEST_EQ(2, m1.rows());
      TEST_EQ(2, m1.cols());
      TEST_EQ(13.0, m1(0, 0));
      TEST_EQ(6.0, m1(1, 0));
      TEST_EQ(69.0, m1(0, 1));
      TEST_EQ(-1.0, m1(1, 1));
      }

    void matrix_multiply_assignment2()
      {
      mat m1(2, 3);
      mat m2(3, 2);
      mat m2zero(3, 2);

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      m1 *= (m2 + m2zero + m2zero);
      TEST_EQ(2, m1.rows());
      TEST_EQ(2, m1.cols());
      TEST_EQ(13.0, m1(0, 0));
      TEST_EQ(6.0, m1(1, 0));
      TEST_EQ(69.0, m1(0, 1));
      TEST_EQ(-1.0, m1(1, 1));
      }

    struct matrix_add_asignment_array : public matrix_fixture {
      void test()
        {
        matf6 res(m16);
        res += m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(2.f, res(0, 1));
        TEST_EQ(4.f, res(0, 2));
        TEST_EQ(6.f, res(1, 0));
        TEST_EQ(8.f, res(1, 1));
        TEST_EQ(10.f, res(1, 2));

        res += (m + m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(4.f, res(0, 1));
        TEST_EQ(8.f, res(0, 2));
        TEST_EQ(12.f, res(1, 0));
        TEST_EQ(16.f, res(1, 1));
        TEST_EQ(20.f, res(1, 2));
        }
      };

    struct matrix_subtract_assignment_array : public matrix_fixture {
      void test()
        {
        matf6 res(2, 3);
        for (auto& v : res)
          v = 0;
        res -= m;
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-1.f, res(0, 1));
        TEST_EQ(-2.f, res(0, 2));
        TEST_EQ(-3.f, res(1, 0));
        TEST_EQ(-4.f, res(1, 1));
        TEST_EQ(-5.f, res(1, 2));
        res -= (m + m);
        TEST_EQ(0.f, res(0, 0));
        TEST_EQ(-3.f, res(0, 1));
        TEST_EQ(-6.f, res(0, 2));
        TEST_EQ(-9.f, res(1, 0));
        TEST_EQ(-12.f, res(1, 1));
        TEST_EQ(-15.f, res(1, 2));
        }
      };

    void matrix_multiply_assignment_array()
      {
      mat6 m1(2, 3);
      mat8 m2(3, 2);

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      m1 *= m2;
      TEST_EQ(2, m1.rows());
      TEST_EQ(2, m1.cols());
      TEST_EQ(13.0, m1(0, 0));
      TEST_EQ(6.0, m1(1, 0));
      TEST_EQ(69.0, m1(0, 1));
      TEST_EQ(-1.0, m1(1, 1));
      }

    void matrix_multiply_assignment2_array()
      {
      mat6 m1(2, 3);
      mat6 m2(3, 2);
      mat6 m2zero(3, 2);
      for (auto& v : m2zero)
        v = 0;

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      m1 *= (m2 + m2zero + m2zero);
      TEST_EQ(2, m1.rows());
      TEST_EQ(2, m1.cols());
      TEST_EQ(13.0, m1(0, 0));
      TEST_EQ(6.0, m1(1, 0));
      TEST_EQ(69.0, m1(0, 1));
      TEST_EQ(-1.0, m1(1, 1));
      }

    struct matrix_diagonal : public matrix_fixture {
      void test()
        {
        matf d = diagonal(m);
        TEST_EQ(2, d.rows());
        TEST_EQ(1, d.cols());
        TEST_EQ(0.f, d(0));
        TEST_EQ(4.f, d(1));

        d = diagonal(m + m16);
        TEST_EQ(2, d.rows());
        TEST_EQ(1, d.cols());
        TEST_EQ(0.f, d(0));
        TEST_EQ(8.f, d(1));

        d = 2.f*diagonal(m + m16) + 1.f;
        TEST_EQ(2, d.rows());
        TEST_EQ(1, d.cols());
        TEST_EQ(1.f, d(0));
        TEST_EQ(17.f, d(1));
        }
      };

    void block_matrix()
      {
      mat m(4, 4);
      m << 1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16;
      mat b = block(m, 1, 1, 2, 2);
      TEST_EQ(6.f, b(0, 0));
      TEST_EQ(7.f, b(0, 1));
      TEST_EQ(10.f, b(1, 0));
      TEST_EQ(11.f, b(1, 1));

      b = block(m, 0, 0, 3, 3);
      TEST_EQ(1.f, b(0, 0));
      TEST_EQ(2.f, b(0, 1));
      TEST_EQ(3.f, b(0, 2));
      TEST_EQ(5, b(1, 0));
      TEST_EQ(6, b(1, 1));
      TEST_EQ(7, b(1, 2));
      TEST_EQ(9, b(2, 0));
      TEST_EQ(10.f, b(2, 1));
      TEST_EQ(11.f, b(2, 2));

      mat d = diagonal(block(m, 0, 0, 3, 3));
      TEST_EQ(1.f, d(0));
      TEST_EQ(6.f, d(1));
      TEST_EQ(11.f, d(2));

      b = block(block(m, 0, 0, 3, 3), 1, 1, 2, 2);
      TEST_EQ(6.f, b(0, 0));
      TEST_EQ(7.f, b(0, 1));
      TEST_EQ(10.f, b(1, 0));
      TEST_EQ(11.f, b(1, 1));
      }

    struct matrix_comparison : public matrix_fixture {
      void test()
        {
        mat m3(3, 2);
        TEST_ASSERT(m == m16);
        TEST_ASSERT(!(m != m16));
        TEST_ASSERT(m != m3);
        TEST_ASSERT(!(m == m3));

        TEST_ASSERT(m == (m16 + m3));
        TEST_ASSERT(!(m != (m16 + m3)));
        TEST_ASSERT(m != (m16 + 1.f));
        TEST_ASSERT(!(m == (m16 + 1.f)));

        TEST_ASSERT((m16 + m3) == m);
        TEST_ASSERT(!((m16 + m3) != m));
        TEST_ASSERT((m16 + 1.f) != m);
        TEST_ASSERT(!((m16 + 1.f) == m));

        TEST_ASSERT(2.f*m16 == m * 2.f);
        TEST_ASSERT(!(2.f*m16 != m * 2.f));
        TEST_ASSERT(2.f*m16 != m * 3.f);
        TEST_ASSERT(!(2.f*m16 == m * 3.f));
        }
      };

    struct matrix_capacity : public matrix_fixture {
      void test()
        {
        TEST_EQ(16, m16.capacity());
        TEST_ASSERT(m.capacity() > 100000);
        }
      };

    void matrix_init_double()
      {
      mat m = zeros(3, 4);
      TEST_EQ(3, m.rows());
      TEST_EQ(4, m.cols());
      for (auto value : m)
        TEST_EQ(0.f, value);

      m = ones(5, 7);
      TEST_EQ(5, m.rows());
      TEST_EQ(7, m.cols());
      for (auto value : m)
        TEST_EQ(1.f, value);

      m = identity(2, 3);
      TEST_EQ(2, m.rows());
      TEST_EQ(3, m.cols());
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, m(i, j));

      mat16 m16 = zeros(3, 4);
      TEST_EQ(3, m16.rows());
      TEST_EQ(4, m16.cols());
      for (auto value : m16)
        TEST_EQ(0.f, value);

      m16 = ones(1, 2);
      TEST_EQ(1, m16.rows());
      TEST_EQ(2, m16.cols());
      for (auto value : m16)
        TEST_EQ(1.f, value);

      m16 = identity(2, 3);
      TEST_EQ(2, m16.rows());
      TEST_EQ(3, m16.cols());
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, m16(i, j));
      }

    void matrix_init_float()
      {
      matf m = zeros<float>(3, 4);
      TEST_EQ(3, m.rows());
      TEST_EQ(4, m.cols());
      for (auto value : m)
        TEST_EQ(0.f, value);

      m = ones<float>(5, 7);
      TEST_EQ(5, m.rows());
      TEST_EQ(7, m.cols());
      for (auto value : m)
        TEST_EQ(1.f, value);

      m = identity<float>(2, 3);
      TEST_EQ(2, m.rows());
      TEST_EQ(3, m.cols());
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, m(i, j));

      matf16 m16 = zeros<float>(3, 4);
      TEST_EQ(3, m16.rows());
      TEST_EQ(4, m16.cols());
      for (auto value : m16)
        TEST_EQ(0.f, value);

      m16 = ones<float>(1, 2);
      TEST_EQ(1, m16.rows());
      TEST_EQ(2, m16.cols());
      for (auto value : m16)
        TEST_EQ(1.f, value);

      m16 = identity<float>(2, 3);
      TEST_EQ(2, m16.rows());
      TEST_EQ(3, m16.cols());
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, m16(i, j));
      }

    void lu_dcmp()
      {
      mat m(2, 2);
      m << 3, 4, 3, 6;
      mat m2 = m;
      std::vector<uint64_t> permutations;
      double d;
      ludcmp(m, permutations, d);
      TEST_EQ(0, permutations[0]);
      TEST_EQ(1, permutations[1]);
      mat l = m;
      mat u = m;
      u(1, 0) = 0.0;
      l(0, 1) = 0.0;
      l(0, 0) = 1.0;
      l(1, 1) = 1.0;
      mat m3 = l * u;
      TEST_EQ(m3(0, 0), m2(0, 0));
      TEST_EQ(m3(0, 1), m2(0, 1));
      TEST_EQ(m3(1, 0), m2(1, 0));
      TEST_EQ(m3(1, 1), m2(1, 1));
      }

    void lu_dcmp2()
      {
      mat m(2, 2);
      m << 4, 3, 6, 3;
      mat m2 = m;
      std::vector<uint64_t> permutations;
      double d;
      ludcmp(m, permutations, d);
      TEST_EQ(1, permutations[0]);
      TEST_EQ(1, permutations[1]);
      mat l = m;
      mat u = m;
      u(1, 0) = 0.0;
      l(0, 1) = 0.0;
      l(0, 0) = 1.0;
      l(1, 1) = 1.0;
      mat m3 = l * u;
      TEST_EQ(m3(0, 0), m2(1, 0));
      TEST_EQ(m3(0, 1), m2(1, 1));
      TEST_EQ(m3(1, 0), m2(0, 0));
      TEST_EQ(m3(1, 1), m2(0, 1));
      }

    void lu_bksb()
      {
      mat m(2, 2);
      m << 4, 3, 6, 3;
      mat v(2, 1);
      v << 10, 12;
      std::vector<uint64_t> permutations;
      double d;
      ludcmp(m, permutations, d);
      lubksb(v, m, permutations);
      TEST_EQ(1.0, v(0));
      TEST_EQ(2.0, v(1));
      }


    void lu_solve()
      {
      mat m(2, 2);
      m << 4, 3, 6, 3;
      mat b(2);
      b << 10, 12;
      solve(m, b);
      TEST_EQ(1, b(0));
      TEST_EQ(2, b(1));
      }

    void lu_invert()
      {
      mat m(2, 2), invm;
      m << 4, 3, 6, 3;
      mat b(2);
      b << 10, 12;
      invert(invm, m);
      mat x = invm * b;
      TEST_EQ(1, x(0));
      TEST_EQ(2, x(1));
      }

    void lu_determinant()
      {
      mat m(2, 2);
      m << 4, 3, 6, 3;
      double d = determinant(m);
      TEST_EQ(-6.0, d);
      }

    void lu_dcmp_array()
      {
      matf16 m(2, 2);
      m << 3, 4, 3, 6;
      matf16 m2 = m;
      std::vector<uint64_t> permutations;
      float d;
      ludcmp(m, permutations, d);
      TEST_EQ(0, permutations[0]);
      TEST_EQ(1, permutations[1]);
      matf16 l = m;
      matf16 u = m;
      u(1, 0) = 0.0;
      l(0, 1) = 0.0;
      l(0, 0) = 1.0;
      l(1, 1) = 1.0;
      matf16 m3 = l * u;
      TEST_EQ(m3(0, 0), m2(0, 0));
      TEST_EQ(m3(0, 1), m2(0, 1));
      TEST_EQ(m3(1, 0), m2(1, 0));
      TEST_EQ(m3(1, 1), m2(1, 1));
      }

    void lu_dcmp2_array()
      {
      matf16 m(2, 2);
      m << 4, 3, 6, 3;
      matf16 m2 = m;
      std::vector<uint64_t> permutations;
      float d;
      ludcmp(m, permutations, d);
      TEST_EQ(1, permutations[0]);
      TEST_EQ(1, permutations[1]);
      matf16 l = m;
      matf16 u = m;
      u(1, 0) = 0.0;
      l(0, 1) = 0.0;
      l(0, 0) = 1.0;
      l(1, 1) = 1.0;
      matf16 m3 = l * u;
      TEST_EQ(m3(0, 0), m2(1, 0));
      TEST_EQ(m3(0, 1), m2(1, 1));
      TEST_EQ(m3(1, 0), m2(0, 0));
      TEST_EQ(m3(1, 1), m2(0, 1));
      }

    void lu_bksb_array()
      {
      matf16 m(2, 2);
      m << 4, 3, 6, 3;
      matf2 v(2, 1);
      v << 10, 12;
      std::vector<uint64_t> permutations;
      float d;
      ludcmp(m, permutations, d);
      lubksb(v, m, permutations);
      TEST_EQ(1.f, v(0));
      TEST_EQ(2.f, v(1));
      }


    void lu_solve_array()
      {
      matf16 m(2, 2);
      m << 4, 3, 6, 3;
      matf2 b(2);
      b << 10, 12;
      solve(m, b);
      TEST_EQ(1, b(0));
      TEST_EQ(2, b(1));
      }

    void lu_invert_array()
      {
      matf16 m(2, 2), invm;
      m << 4, 3, 6, 3;
      matf16 b(2);
      b << 10, 12;
      invert(invm, m);
      matf16 x = invm * b;
      TEST_EQ(1, x(0));
      TEST_EQ(2, x(1));
      }

    void lu_determinant_array()
      {
      matf16 m(2, 2);
      m << 4, 3, 6, 3;
      float d = determinant(m);
      TEST_EQ(-6.f, d);
      }

    void cholesky1()
      {
      mat m(3, 3);
      m << 2, 2, 4, 2, 4, 6, 4, 6, 18;
      mat old_m = m;

      cholesky(m);

      mat u(m.rows(), m.cols());
      for (size_t i = 0; i < u.rows(); ++i)
        for (size_t j = i; j < u.cols(); ++j)
          u(i, j) = m(i, j);
      mat m2 = transpose(u)*u;
      for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
          TEST_EQ_CLOSE(old_m(i, j), m2(i, j), 0.00000001);
      }

    void cholesky2()
      {
      mat m = identity(5, 5)*2.0;
      for (uint64_t i = 0; i < 5; ++i)
        {
        if (i > 0)
          m(i, i - 1) = -1.0;
        if (i < 4)
          m(i, i + 1) = -1.0;
        }

      mat old_m = m;

      cholesky(m);

      mat u(m.rows(), m.cols());
      for (size_t i = 0; i < u.rows(); ++i)
        for (size_t j = i; j < u.cols(); ++j)
          u(i, j) = m(i, j);
      mat m2 = transpose(u)*u;
      for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
          TEST_EQ_CLOSE(old_m(i, j), m2(i, j), 0.00000001);
      }

    void back_substitution()
      {
      mat u(3, 3);
      u(0, 0) = 1.0;
      u(0, 1) = 2.0;
      u(0, 2) = 3.0;
      u(1, 1) = 4.0;
      u(1, 2) = 5.0;
      u(2, 2) = 6.0;
      mat x(3);
      mat b(3);
      b(0) = 1.0;
      b(1) = 2.0;
      b(2) = 3.0;
      back_substitution(x, b, u);
      mat c = u * x;
      for (uint64_t i = 0; i < b.rows(); ++i)
        TEST_EQ_CLOSE(b(i), c(i), 0.00000001);
      }

    void forward_substitution()
      {
      mat l(3, 3);
      l(0, 0) = 1.0;
      l(1, 0) = 2.0;
      l(2, 0) = 3.0;
      l(1, 1) = 4.0;
      l(2, 1) = 5.0;
      l(2, 2) = 6.0;
      mat x(3);
      mat b(3);
      b(0) = 1.0;
      b(1) = 2.0;
      b(2) = 3.0;
      forward_substitution(x, b, l);
      mat c = l * x;
      for (uint64_t i = 0; i < b.rows(); ++i)
        TEST_EQ_CLOSE(b(i), c(i), 0.00000001);
      }

    void solve_cholesky()
      {
      mat m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      mat b(3, 1);
      b << -20, -43, 192;
      cholesky(m);
      mat x, y;
      forward_substitution(y, b, mat(transpose(m)));
      back_substitution(x, y, m);
      TEST_EQ(1.0, x(0));
      TEST_EQ(2.0, x(1));
      TEST_EQ(3.0, x(2));
      }

    void solve_cholesky_2()
      {
      mat m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      mat b(3, 1);
      b << -20, -43, 192;
      solve_cholesky(m, b);
      TEST_EQ(1.0, b(0));
      TEST_EQ(2.0, b(1));
      TEST_EQ(3.0, b(2));
      }

    void cholesky1_array()
      {
      matf16 m(3, 3);
      m << 2, 2, 4, 2, 4, 6, 4, 6, 18;
      matf16 old_m = m;

      cholesky(m);

      matf16 u = zeros(m.rows(), m.cols());
      for (size_t i = 0; i < u.rows(); ++i)
        for (size_t j = i; j < u.cols(); ++j)
          u(i, j) = m(i, j);
      matf16 m2 = transpose(u)*u;
      for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
          TEST_EQ_CLOSE(old_m(i, j), m2(i, j), 0.00001f);
      }

    void cholesky2_array()
      {
      matf25 m = identity(5, 5)*2.f;
      for (uint64_t i = 0; i < 5; ++i)
        {
        if (i > 0)
          m(i, i - 1) = -1.0;
        if (i < 4)
          m(i, i + 1) = -1.0;
        }

      matf25 old_m = m;

      cholesky(m);

      matf25 u = zeros(m.rows(), m.cols());
      for (size_t i = 0; i < u.rows(); ++i)
        for (size_t j = i; j < u.cols(); ++j)
          u(i, j) = m(i, j);
      matf25 m2 = transpose(u)*u;
      for (size_t i = 0; i < m.rows(); ++i)
        for (size_t j = 0; j < m.cols(); ++j)
          TEST_EQ_CLOSE(old_m(i, j), m2(i, j), 0.00001f);
      }

    void back_substitution_array()
      {
      matf9 u = zeros(3, 3);
      u(0, 0) = 1.0;
      u(0, 1) = 2.0;
      u(0, 2) = 3.0;
      u(1, 1) = 4.0;
      u(1, 2) = 5.0;
      u(2, 2) = 6.0;
      matf3 x(3);
      matf3 b(3);
      b(0) = 1.0;
      b(1) = 2.0;
      b(2) = 3.0;
      back_substitution(x, b, u);
      matf3 c = u * x;
      for (uint64_t i = 0; i < b.rows(); ++i)
        TEST_EQ_CLOSE(b(i), c(i), 0.00001f);
      }

    void forward_substitution_array()
      {
      matf9 l = zeros(3, 3);
      l(0, 0) = 1.0;
      l(1, 0) = 2.0;
      l(2, 0) = 3.0;
      l(1, 1) = 4.0;
      l(2, 1) = 5.0;
      l(2, 2) = 6.0;
      matf3 x(3);
      matf3 b(3);
      b(0) = 1.0;
      b(1) = 2.0;
      b(2) = 3.0;
      forward_substitution(x, b, l);
      matf3 c = l * x;
      for (uint64_t i = 0; i < b.rows(); ++i)
        TEST_EQ_CLOSE(b(i), c(i), 0.00001f);
      }

    void solve_cholesky_array()
      {
      matf9 m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      matf3 b(3, 1);
      b << -20, -43, 192;
      cholesky(m);
      matf3 x, y;
      forward_substitution(y, b, matf16(transpose(m)));
      back_substitution(x, y, m);
      TEST_EQ(1.f, x(0));
      TEST_EQ(2.f, x(1));
      TEST_EQ(3.f, x(2));
      }

    void solve_cholesky_2_array()
      {
      matf9 m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      matf3 b(3, 1);
      b << -20, -43, 192;
      solve_cholesky(m, b);
      TEST_EQ(1.f, b(0));
      TEST_EQ(2.f, b(1));
      TEST_EQ(3.f, b(2));
      }

    void very_large_size_expression_uncomputed()
      {
      auto matrix_expr = identity(0xffffffffffffffff, 0xffffffffffffffff)*3.0 - 1.0;
      TEST_EQ(0xffffffffffffffff, matrix_expr.rows());
      TEST_EQ(0xffffffffffffffff, matrix_expr.cols());
      }

    void solve_qr()
      {
      mat m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      mat b(3, 1);
      b << -20, -43, 192;
      solve_qr(m, b);
      TEST_EQ_CLOSE(1.0, b(0), 1e-12);
      TEST_EQ_CLOSE(2.0, b(1), 1e-12);
      TEST_EQ_CLOSE(3.0, b(2), 1e-12);
      }

    void qr()
      {
      mat m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      mat m2 = m;
      mat qt;
      std::vector<uint64_t> p;
      qrdcmp(m, qt, p, false);
      mat res = transpose(qt)*m;
      for (int i = 0; i < res.rows(); ++i)
        for (int j = 0; j < res.cols(); ++j)
          TEST_EQ_CLOSE(m2(i, j), res(i, j), 1e-12);
      mat qqt = transpose(qt)*qt;
      for (int i = 0; i < qqt.rows(); ++i)
        for (int j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-12);
      }

    void qr_rectangular()
      {
      mat m(4, 3);
      m << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      mat m2 = m;
      mat qt;
      std::vector<uint64_t> p;
      qrdcmp(m, qt, p, false);
      mat res = transpose(qt)*m;
      for (int i = 0; i < res.rows(); ++i)
        for (int j = 0; j < res.cols(); ++j)
          TEST_EQ_CLOSE(m2(i, j), res(i, j), 1e-12);
      mat qqt = transpose(qt)*qt;
      for (int i = 0; i < qqt.rows(); ++i)
        for (int j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-12);
      }

    void qr_rectangular2()
      {
      mat m(3, 4);
      m << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      mat m2 = m;
      mat qt;
      std::vector<uint64_t> p;
      qrdcmp(m, qt, p, false);
      mat res = transpose(qt)*m;
      for (int i = 0; i < res.rows(); ++i)
        for (int j = 0; j < res.cols(); ++j)
          TEST_EQ_CLOSE(m2(i, j), res(i, j), 1e-12);
      mat qqt = transpose(qt)*qt;
      for (int i = 0; i < qqt.rows(); ++i)
        for (int j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-12);
      }

    void solve_qr_test_least_squares()
      {
      mat A(3, 2);
      A << 2, 0, -1, 1, 0, 2;
      mat b(3, 1);
      b << 1, 0, -1;
      solve_qr(A, b);
      TEST_EQ(2, b.rows());
      TEST_EQ(1, b.cols());
      TEST_EQ_CLOSE(1.0 / 3.0, b(0), 1e-12);
      TEST_EQ_CLOSE(-1.0 / 3.0, b(1), 1e-12);
      }

    void solve_qr_test_least_squares_2()
      {
      mat A(4, 3);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = double(i);
      mat b(4);
      for (int i = 0; i < 4; ++i)
        b(i) = double(3 * i + 2);
      solve_qr(A, b);
      TEST_EQ_CLOSE(0.0, b(0), 1e-12);
      TEST_EQ_CLOSE(1.0, b(1), 1e-12);
      TEST_EQ_CLOSE(0.0, b(2), 1e-12);
      }

    void solve_qr_array()
      {
      matf9 m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      matf3 b(3, 1);
      b << -20, -43, 192;
      solve_qr(m, b);
      TEST_EQ_CLOSE(1.f, b(0), 1e-2);
      TEST_EQ_CLOSE(2.f, b(1), 1e-2);
      TEST_EQ_CLOSE(3.f, b(2), 1e-2);
      }

    void qr_array()
      {
      matf9 m(3, 3);
      m << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      matf9 m2 = m;
      matf9 qt;
      std::vector<uint64_t> p;
      qrdcmp(m, qt, p, false);
      matf16 res = transpose(qt)*m;
      for (int i = 0; i < res.rows(); ++i)
        for (int j = 0; j < res.cols(); ++j)
          TEST_EQ_CLOSE(m2(i, j), res(i, j), 1e-5);
      matf16 qqt = transpose(qt)*qt;
      for (int i = 0; i < qqt.rows(); ++i)
        for (int j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-5);
      }

    void qr_rectangular_array()
      {
      matf16 m(4, 3);
      m << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      matf16 m2 = m;
      matf16 qt;
      std::vector<uint64_t> p;
      qrdcmp(m, qt, p, false);
      matf16 res = transpose(qt)*m;
      for (int i = 0; i < res.rows(); ++i)
        for (int j = 0; j < res.cols(); ++j)
          TEST_EQ_CLOSE(m2(i, j), res(i, j), 1e-5);
      matf16 qqt = transpose(qt)*qt;
      for (int i = 0; i < qqt.rows(); ++i)
        for (int j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-5);
      }

    void qr_rectangular2_array()
      {
      matf16 m(3, 4);
      m << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      matf16 m2 = m;
      matf16 qt;
      std::vector<uint64_t> p;
      qrdcmp(m, qt, p, false);
      matf16 res = transpose(qt)*m;
      for (int i = 0; i < res.rows(); ++i)
        for (int j = 0; j < res.cols(); ++j)
          TEST_EQ_CLOSE(m2(i, j), res(i, j), 1e-5);
      matf16 qqt = transpose(qt)*qt;
      for (int i = 0; i < qqt.rows(); ++i)
        for (int j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-5);
      }

    void solve_qr_test_least_squares_array()
      {
      matf6 A(3, 2);
      A << 2, 0, -1, 1, 0, 2;
      matf3 b(3, 1);
      b << 1, 0, -1;
      solve_qr(A, b);
      TEST_EQ(2, b.rows());
      TEST_EQ(1, b.cols());
      TEST_EQ_CLOSE(1.f / 3.f, b(0), 1e-5);
      TEST_EQ_CLOSE(-1.f / 3.f, b(1), 1e-5);
      }

    void solve_qr_test_least_squares_2_array()
      {
      matf12 A(4, 3);
      for (int i = 1; i <= 12; ++i)
        A[(i - 1) / 3][(i - 1) % 3] = float(i);
      matf4 b(4);
      for (int i = 0; i < 4; ++i)
        b(i) = float(3 * i + 2);
      solve_qr(A, b);
      TEST_EQ_CLOSE(0.f, b(0), 1e-5);
      TEST_EQ_CLOSE(1.f, b(1), 1e-5);
      TEST_EQ_CLOSE(0.f, b(2), 1e-5);
      }

    void norm_tests()
      {
      mat v(5, 1);
      v << 1, 2, 3, 4, 5;
      TEST_EQ(std::sqrt(1 + 4 + 9 + 16 + 25), norm(v));
      mat m(2, 2);
      m << 2, 3, 7, 8;
      TEST_EQ(std::sqrt(4 + 9 + 49 + 64), norm(m));
      mat m2(2, 2);
      m2 << -1, -2, -6, -7;
      TEST_EQ(2.0, norm(m + m2));
      }

    void blocknorm()
      {
      mat m(3, 4);
      m << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      double n = norm(block(m, 0, 1, 3, 1));
      TEST_EQ(std::sqrt(1.0 + 4.0 + 1.0), n);
      }

    void norm_tests_array()
      {
      matf16 v(5, 1);
      v << 1, 2, 3, 4, 5;
      TEST_EQ_CLOSE(std::sqrt(1 + 4 + 9 + 16 + 25), norm(v), 1e-5);
      matf16 m(2, 2);
      m << 2, 3, 7, 8;
      TEST_EQ_CLOSE(std::sqrt(4 + 9 + 49 + 64), norm(m), 1e-5);
      matf16 m2(2, 2);
      m2 << -1, -2, -6, -7;
      TEST_EQ(2.0, norm(m + m2));
      }

    void blocknorm_array()
      {
      matf16 m(3, 4);
      m << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      double n = norm(block(m, 0, 1, 3, 1));
      TEST_EQ(std::sqrt(1.0 + 4.0 + 1.0), n);
      }

    void qrfac()
      {
      uint64_t m = 4;
      uint64_t n = 3;
      mat a(m, n);
      a << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      mat a2 = a;
      std::vector<uint64_t> permutations;
      mat rdiag;
      mat acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);
      mat P = zeros(n, n);
      for (uint64_t j = 0; j < n; ++j)
        P(permutations[j], j) = 1.0;
      mat qt(m, m);
      mat r = a;
      uint64_t mn = m > n ? n : m;
      uint64_t i, j, k;
      for (i = 0; i < m; ++i)
        {
        for (j = 0; j < m; ++j)
          qt(i, j) = 0.0;
        qt(i, i) = 1.0;
        }

      for (k = 0; k < mn; ++k)
        {
        if (a(k, k) != 0.0)
          {
          for (j = 0; j < m; ++j)
            {
            double sum = 0.0;
            for (i = k; i < m; ++i)
              sum += a(i, k) * qt(i, j);
            sum /= a(k, k);
            for (i = k; i < m; ++i)
              qt(i, j) -= sum * a(i, k);
            }
          }
        }

      for (i = 0; i < mn; ++i)
        {
        r(i, i) = rdiag(i);
        for (j = 0; j < i; ++j)
          r(i, j) = 0.0;
        }
      for (i = mn; i < m; ++i)
        {
        for (j = 0; j < n; ++j)
          r(i, j) = 0.0;
        }

      mat right = transpose(qt)*r;
      mat left = a2 * P;
      for (i = 0; i < right.rows(); ++i)
        for (j = 0; j < right.cols(); ++j)
          TEST_EQ_CLOSE(left(i, j), right(i, j), 1e-12);
      mat qqt = transpose(qt)*qt;
      for (i = 0; i < qqt.rows(); ++i)
        for (j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-12);

      }

    void qrfac2()
      {
      uint64_t m = 3;
      uint64_t n = 4;
      mat a(m, n);
      a << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      mat a2 = a;
      std::vector<uint64_t> permutations;
      mat rdiag;
      mat acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);
      mat P = zeros(n, n);
      for (uint64_t j = 0; j < n; ++j)
        P(permutations[j], j) = 1.0;
      mat qt(m, m);
      mat r = a;
      uint64_t mn = m > n ? n : m;
      uint64_t i, j, k;
      for (i = 0; i < m; ++i)
        {
        for (j = 0; j < m; ++j)
          qt(i, j) = 0.0;
        qt(i, i) = 1.0;
        }

      for (k = 0; k < mn; ++k)
        {
        if (a(k, k) != 0.0)
          {
          for (j = 0; j < m; ++j)
            {
            double sum = 0.0;
            for (i = k; i < m; ++i)
              sum += a(i, k) * qt(i, j);
            sum /= a(k, k);
            for (i = k; i < m; ++i)
              qt(i, j) -= sum * a(i, k);
            }
          }
        }

      for (i = 0; i < mn; ++i)
        {
        r(i, i) = rdiag(i);
        for (j = 0; j < i; ++j)
          r(i, j) = 0.0;
        }
      for (i = mn; i < m; ++i)
        {
        for (j = 0; j < n; ++j)
          r(i, j) = 0.0;
        }

      mat right = transpose(qt)*r;
      mat left = a2 * P;
      for (i = 0; i < right.rows(); ++i)
        for (j = 0; j < right.cols(); ++j)
          TEST_EQ_CLOSE(left(i, j), right(i, j), 1e-12);
      mat qqt = transpose(qt)*qt;
      for (i = 0; i < qqt.rows(); ++i)
        for (j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.0 : 0.0, qqt(i, j), 1e-12);

      }


    void qrsolvtest()
      {
      int m = 3;
      int n = 3;
      mat a(3, 3);
      a << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      mat b(3, 1);
      b << -20, -43, 192;

      std::vector<uint64_t> permutations;
      mat rdiag;
      mat acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);

      mat r = a;
      for (int i = 0; i < n; ++i)
        {
        r(i, i) = rdiag(i);
        for (int j = 0; j < i; ++j)
          r(i, j) = 0.0;
        }
      int mn = m > n ? n : m;
      mat qt = identity(m, m);
      for (int k = 0; k < mn; ++k)
        {
        if (a(k, k) != 0.0)
          {
          for (int j = 0; j < m; ++j)
            {
            double sum = 0.0;
            for (int i = k; i < m; ++i)
              sum += a(i, k) * qt(i, j);
            sum /= a(k, k);
            for (int i = k; i < m; ++i)
              qt(i, j) -= sum * a(i, k);
            }
          }
        }

      mat qtbtest = qt * b;

      mat qtb(3, 1);

      for (int j = 0; j < n; ++j)
        {
        if (a[j][j] != 0.0)
          {
          double sum = 0.0;
          for (int i = j; i < m; ++i)
            sum += a[i][j] * b(i);
          double temp = -sum / a[j][j];
          for (int i = j; i < m; ++i)
            b(i) += a[i][j] * temp;
          }
        qtb(j) = b(j);
        }

      TEST_ASSERT(norm(qtbtest - qtb) < 1e-12);

      mat x, sdiag;
      mat diag = zeros(n, 1);

      qrsolv(r, permutations, diag, qtb, x, sdiag);

      TEST_EQ_CLOSE(1.0, x(0), 1e-12);
      TEST_EQ_CLOSE(2.0, x(1), 1e-12);
      TEST_EQ_CLOSE(3.0, x(2), 1e-12);
      }

    void qrsolvtest2()
      {
      int m = 4;
      int n = 3;
      mat a(4, 3);
      for (int i = 1; i <= 12; ++i)
        a[(i - 1) / 3][(i - 1) % 3] = double(i);
      mat b(4);
      for (int i = 0; i < 4; ++i)
        b(i) = double(3 * i + 2);

      mat a2(a);
      mat b2(b);

      std::vector<uint64_t> permutations;
      mat rdiag;
      mat acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);

      mat r = zeros(n, n);
      for (int i = 0; i < n; ++i)
        {
        r(i, i) = rdiag(i);
        for (int j = i + 1; j < n; ++j)
          r(i, j) = a(i, j);
        }
      mat qtb(n, 1);

      for (int j = 0; j < n; ++j)
        {
        if (a[j][j] != 0.0)
          {
          double sum = 0.0;
          for (int i = j; i < m; ++i)
            sum += a[i][j] * b(i);
          double temp = -sum / a[j][j];
          for (int i = j; i < m; ++i)
            b(i) += a[i][j] * temp;
          }
        qtb(j) = b(j);
        }

      mat x, sdiag;
      mat diag = zeros(n, 1);

      qrsolv(r, permutations, diag, qtb, x, sdiag);

      double residu = norm(a2*x - b2);

      TEST_EQ_CLOSE(0.0, residu, 1e-12);
      }


    void qrsolvtest3()
      {
      int m = 3;
      int n = 2;
      mat a(3, 2);
      a << 2, 0, -1, 1, 0, 2;
      mat b(3, 1);
      b << 1, 0, -1;

      mat a2(a);
      mat b2(b);

      std::vector<uint64_t> permutations;
      mat rdiag;
      mat acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);

      mat r = zeros(n, n);
      for (int i = 0; i < n; ++i)
        {
        r(i, i) = rdiag(i);
        for (int j = i + 1; j < n; ++j)
          r(i, j) = a(i, j);
        }
      mat qtb(n, 1);

      for (int j = 0; j < n; ++j)
        {
        if (a[j][j] != 0.0)
          {
          double sum = 0.0;
          for (int i = j; i < m; ++i)
            sum += a[i][j] * b(i);
          double temp = -sum / a[j][j];
          for (int i = j; i < m; ++i)
            b(i) += a[i][j] * temp;
          }
        qtb(j) = b(j);
        }

      mat x, sdiag;
      mat diag = zeros(n, 1);

      qrsolv(r, permutations, diag, qtb, x, sdiag);

      TEST_EQ_CLOSE(1.0 / 3.0, x(0), 1e-12);
      TEST_EQ_CLOSE(-1.0 / 3.0, x(1), 1e-12);
      }

    void qrfac_array()
      {
      uint64_t m = 4;
      uint64_t n = 3;
      matf16 a(m, n);
      a << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      matf16 a2 = a;
      std::vector<uint64_t> permutations;
      matf3 rdiag;
      matf3 acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);
      matf16 P = zeros(n, n);
      for (uint64_t j = 0; j < n; ++j)
        P(permutations[j], j) = 1.0;
      matf16 qt(m, m);
      matf16 r = a;
      uint64_t mn = m > n ? n : m;
      uint64_t i, j, k;
      for (i = 0; i < m; ++i)
        {
        for (j = 0; j < m; ++j)
          qt(i, j) = 0.0;
        qt(i, i) = 1.0;
        }

      for (k = 0; k < mn; ++k)
        {
        if (a(k, k) != 0.f)
          {
          for (j = 0; j < m; ++j)
            {
            float sum = 0.f;
            for (i = k; i < m; ++i)
              sum += a(i, k) * qt(i, j);
            sum /= a(k, k);
            for (i = k; i < m; ++i)
              qt(i, j) -= sum * a(i, k);
            }
          }
        }

      for (i = 0; i < mn; ++i)
        {
        r(i, i) = rdiag(i);
        for (j = 0; j < i; ++j)
          r(i, j) = 0.f;
        }
      for (i = mn; i < m; ++i)
        {
        for (j = 0; j < n; ++j)
          r(i, j) = 0.f;
        }

      matf16 right = transpose(qt)*r;
      matf16 left = a2 * P;
      for (i = 0; i < right.rows(); ++i)
        for (j = 0; j < right.cols(); ++j)
          TEST_EQ_CLOSE(left(i, j), right(i, j), 1e-5);
      matf16 qqt = transpose(qt)*qt;
      for (i = 0; i < qqt.rows(); ++i)
        for (j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.f : 0.f, qqt(i, j), 1e-5);

      }

    void qrfac2_array()
      {
      uint64_t m = 3;
      uint64_t n = 4;
      matf16 a(m, n);
      a << 1, -1, 0, 1, -1, 2, 1, 1, 0, 1, 1, 0;
      matf16 a2 = a;
      std::vector<uint64_t> permutations;
      matf4 rdiag;
      matf4 acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);
      matf16 P = zeros(n, n);
      for (uint64_t j = 0; j < n; ++j)
        P(permutations[j], j) = 1.f;
      matf16 qt(m, m);
      matf16 r = a;
      uint64_t mn = m > n ? n : m;
      uint64_t i, j, k;
      for (i = 0; i < m; ++i)
        {
        for (j = 0; j < m; ++j)
          qt(i, j) = 0.f;
        qt(i, i) = 1.f;
        }

      for (k = 0; k < mn; ++k)
        {
        if (a(k, k) != 0.f)
          {
          for (j = 0; j < m; ++j)
            {
            float sum = 0.f;
            for (i = k; i < m; ++i)
              sum += a(i, k) * qt(i, j);
            sum /= a(k, k);
            for (i = k; i < m; ++i)
              qt(i, j) -= sum * a(i, k);
            }
          }
        }

      for (i = 0; i < mn; ++i)
        {
        r(i, i) = rdiag(i);
        for (j = 0; j < i; ++j)
          r(i, j) = 0.f;
        }
      for (i = mn; i < m; ++i)
        {
        for (j = 0; j < n; ++j)
          r(i, j) = 0.f;
        }

      matf16 right = transpose(qt)*r;
      matf16 left = a2 * P;
      for (i = 0; i < right.rows(); ++i)
        for (j = 0; j < right.cols(); ++j)
          TEST_EQ_CLOSE(left(i, j), right(i, j), 1e-5);
      mat qqt = transpose(qt)*qt;
      for (i = 0; i < qqt.rows(); ++i)
        for (j = 0; j < qqt.cols(); ++j)
          TEST_EQ_CLOSE(i == j ? 1.f : 0.f, qqt(i, j), 1e-5);

      }

    void qrsolvtest_array()
      {
      int m = 3;
      int n = 3;
      matf16 a(3, 3);
      a << 4, 12, -16, 12, 37, -43, -16, -43, 98;
      matf3 b(3, 1);
      b << -20, -43, 192;

      std::vector<uint64_t> permutations;
      matf3 rdiag;
      matf3 acnorm;
      qrfac(a, true, permutations, rdiag, acnorm);

      matf16 r = a;
      for (int i = 0; i < n; ++i)
        {
        r(i, i) = rdiag(i);
        for (int j = 0; j < i; ++j)
          r(i, j) = 0.f;
        }
      int mn = m > n ? n : m;
      matf16 qt = identity(m, m);
      for (int k = 0; k < mn; ++k)
        {
        if (a(k, k) != 0.f)
          {
          for (int j = 0; j < m; ++j)
            {
            float sum = 0.f;
            for (int i = k; i < m; ++i)
              sum += a(i, k) * qt(i, j);
            sum /= a(k, k);
            for (int i = k; i < m; ++i)
              qt(i, j) -= sum * a(i, k);
            }
          }
        }

      matf3 qtbtest = qt * b;

      matf3 qtb(3, 1);

      for (int j = 0; j < n; ++j)
        {
        if (a[j][j] != 0.f)
          {
          float sum = 0.f;
          for (int i = j; i < m; ++i)
            sum += a[i][j] * b(i);
          float temp = -sum / a[j][j];
          for (int i = j; i < m; ++i)
            b(i) += a[i][j] * temp;
          }
        qtb(j) = b(j);
        }

      TEST_ASSERT(norm(qtbtest - qtb) < 1e-5);

      matf3 x, sdiag;
      matf3 diag = zeros(n, 1);

      qrsolv(r, permutations, diag, qtb, x, sdiag);

      TEST_EQ_CLOSE(1.f, x(0), 1e-2);
      TEST_EQ_CLOSE(2.f, x(1), 1e-2);
      TEST_EQ_CLOSE(3.f, x(2), 1e-2);
      }


    bool fdjac_test_fun(const mat& x, mat& f, void*)
      {
      f(0) = x(0);
      f(1) = x(1);
      f(2) = x(0)*x(0) + x(1)*x(1);
      return true;
      }


    void fdjac_test()
      {
      mat x = zeros(2, 1);
      mat fvec(3);
      fdjac_test_fun(x, fvec, nullptr);
      mat fjac;
      TEST_ASSERT(fdjac(&fdjac_test_fun, x, fvec, fjac, 1.0));
      TEST_EQ(1.0, fjac(0, 0));
      TEST_EQ(0.0, fjac(0, 1));
      TEST_EQ(0.0, fjac(1, 0));
      TEST_EQ(1.0, fjac(1, 1));
      TEST_EQ(1.0, fjac(2, 0));
      TEST_EQ(1.0, fjac(2, 1));
      TEST_ASSERT(fdjac(&fdjac_test_fun, x, fvec, fjac, 0.5));
      TEST_EQ(1.0, fjac(0, 0));
      TEST_EQ(0.0, fjac(0, 1));
      TEST_EQ(0.0, fjac(1, 0));
      TEST_EQ(1.0, fjac(1, 1));
      TEST_EQ_CLOSE(std::sqrt(0.5), fjac(2, 0), 1e-12);
      TEST_EQ_CLOSE(std::sqrt(0.5), fjac(2, 1), 1e-12);
      }



    bool fdjac_test_array_fun(const matf2& x, matf3& f, void*)
      {
      f(0) = x(0);
      f(1) = x(1);
      f(2) = x(0)*x(0) + x(1)*x(1);
      return true;
      }


    void fdjac_test_array()
      {
      matf2 x = zeros(2, 1);
      matf3 fvec(3);
      fdjac_test_array_fun(x, fvec, nullptr);
      matf6 fjac;
      TEST_ASSERT(fdjac(&fdjac_test_array_fun, x, fvec, fjac, 1.f));
      TEST_EQ(1.0, fjac(0, 0));
      TEST_EQ(0.0, fjac(0, 1));
      TEST_EQ(0.0, fjac(1, 0));
      TEST_EQ(1.0, fjac(1, 1));
      TEST_EQ(1.0, fjac(2, 0));
      TEST_EQ(1.0, fjac(2, 1));
      TEST_ASSERT(fdjac(&fdjac_test_array_fun, x, fvec, fjac, 0.5f));
      TEST_EQ(1.0, fjac(0, 0));
      TEST_EQ(0.0, fjac(0, 1));
      TEST_EQ(0.0, fjac(1, 0));
      TEST_EQ(1.0, fjac(1, 1));
      TEST_EQ_CLOSE(std::sqrt(0.5f), fjac(2, 0), 1e-5);
      TEST_EQ_CLOSE(std::sqrt(0.5f), fjac(2, 1), 1e-5);
      }


    bool error_functional(const mat& solution, mat& f, void* data)
      {
      // y(x) = a * cos (b*x) + b * sin(a * x)    
      int m = 20;
      f.resize(m, 1);
      for (int i = 0; i < m; ++i)
        {
        double a = *((double*)data);
        double b = *(((double*)data) + 1);
        double x = i * 1.0 / (double)m;
        double actual_y = a * std::cos(b*x) + b * std::sin(a * x);
        a = solution(0);
        b = solution(1);
        double current_y = a * std::cos(b*x) + b * std::sin(a * x);
        double error = (actual_y - current_y)*(actual_y - current_y);
        f(i) = error;
        }
      return true;
      }


    void lm_test()
      {
      std::vector<double> coeff;
      coeff.push_back(2.0);
      coeff.push_back(3.0);
      uint64_t info, nfev;
      mat solution, fvec;
      solution = ones(2, 1);
      lmdif0(error_functional, solution, fvec, info, nfev, (void*)coeff.data());
      TEST_EQ_CLOSE(2.0, solution(0), 1e-8);
      TEST_EQ_CLOSE(3.0, solution(1), 1e-8);
      TEST_EQ(1, info);
      }



    bool error_functional_array(const matf2& solution, matf20& f, void* data)
      {
      // y(x) = a * cos (b*x) + b * sin(a * x)    
      int m = 20;
      f.resize(m, 1);
      for (int i = 0; i < m; ++i)
        {
        float a = *((float*)data);
        float b = *(((float*)data) + 1);
        float x = i * 1.f / (float)m;
        float actual_y = a * std::cos(b*x) + b * std::sin(a * x);
        a = solution(0);
        b = solution(1);
        float current_y = a * std::cos(b*x) + b * std::sin(a * x);
        float error = (actual_y - current_y)*(actual_y - current_y);
        f(i) = error;
        }
      return true;

      }

    void lm_test_array()
      {
      std::vector<float> coeff;
      coeff.push_back(2.0);
      coeff.push_back(3.0);
      uint64_t info, nfev;
      matf2 solution;
      matf20 fvec;
      solution = ones(2, 1);
      lmdif0(error_functional_array, solution, fvec, info, nfev, (void*)coeff.data());
      TEST_EQ_CLOSE(2.f, solution(0), 1e-2);
      TEST_EQ_CLOSE(3.f, solution(1), 1e-2);
      TEST_EQ(1, info);
      }

    void tred2test()
      {
      mat matrix(4, 4), orig(4, 4);
      matrix(0, 0) = 4;
      matrix(0, 1) = 1;
      matrix(0, 2) = -2;
      matrix(0, 3) = 2;
      matrix(1, 0) = 1;
      matrix(1, 1) = 2;
      matrix(1, 2) = 0;
      matrix(1, 3) = 1;
      matrix(2, 0) = -2;
      matrix(2, 1) = 0;
      matrix(2, 2) = 3;
      matrix(2, 3) = -2;
      matrix(3, 0) = 2;
      matrix(3, 1) = 1;
      matrix(3, 2) = -2;
      matrix(3, 3) = -1;
      orig = matrix;

      mat diag, subdiag;
      tred2(diag, subdiag, matrix);

      mat result = transpose(matrix)*orig*(matrix);

      TEST_EQ_CLOSE(result(0, 0), diag(0), 1e-8);
      TEST_EQ_CLOSE(result(1, 1), diag(1), 1e-8);
      TEST_EQ_CLOSE(result(2, 2), diag(2), 1e-8);
      TEST_EQ_CLOSE(result(3, 3), diag(3), 1e-8);

      TEST_EQ_CLOSE(0, subdiag(0), 0.0001);
      TEST_EQ_CLOSE(result(1, 0), subdiag(1), 1e-8);
      TEST_EQ_CLOSE(result(2, 1), subdiag(2), 1e-8);
      TEST_EQ_CLOSE(result(3, 2), subdiag(3), 1e-8);

      TEST_EQ_CLOSE(result(0, 1), subdiag(1), 1e-8);
      TEST_EQ_CLOSE(result(1, 2), subdiag(2), 1e-8);
      TEST_EQ_CLOSE(result(2, 3), subdiag(3), 1e-8);

      TEST_EQ_CLOSE(result(0, 2), 0, 1e-8);
      TEST_EQ_CLOSE(result(0, 3), 0, 1e-8);
      TEST_EQ_CLOSE(result(1, 3), 0, 1e-8);
      TEST_EQ_CLOSE(result(2, 0), 0, 1e-8);
      TEST_EQ_CLOSE(result(3, 0), 0, 1e-8);
      TEST_EQ_CLOSE(result(3, 1), 0, 1e-8);
      }


    void tqlitest()
      {
      mat matrix(4, 4), orig(4, 4);
      matrix(0, 0) = 4;
      matrix(0, 1) = 1;
      matrix(0, 2) = -2;
      matrix(0, 3) = 2;
      matrix(1, 0) = 1;
      matrix(1, 1) = 2;
      matrix(1, 2) = 0;
      matrix(1, 3) = 1;
      matrix(2, 0) = -2;
      matrix(2, 1) = 0;
      matrix(2, 2) = 3;
      matrix(2, 3) = -2;
      matrix(3, 0) = 2;
      matrix(3, 1) = 1;
      matrix(3, 2) = -2;
      matrix(3, 3) = -1;
      orig = matrix;

      mat diag, subdiag;
      tred2(diag, subdiag, matrix);
      TEST_ASSERT(tqli(diag, matrix, subdiag));

      mat eigv(4);
      eigv(0) = matrix(0, 0);
      eigv(1) = matrix(1, 0);
      eigv(2) = matrix(2, 0);
      eigv(3) = matrix(3, 0);

      mat result;
      result = orig * eigv;
      double eig = diag(0);

      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);

      eigv(0) = matrix(0, 1);
      eigv(1) = matrix(1, 1);
      eigv(2) = matrix(2, 1);
      eigv(3) = matrix(3, 1);
      result = orig * eigv;
      eig = diag(1);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);

      eigv(0) = matrix(0, 2);
      eigv(1) = matrix(1, 2);
      eigv(2) = matrix(2, 2);
      eigv(3) = matrix(3, 2);
      result = orig * eigv;
      eig = diag(2);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);

      eigv(0) = matrix(0, 3);
      eigv(1) = matrix(1, 3);
      eigv(2) = matrix(2, 3);
      eigv(3) = matrix(3, 3);
      result = orig * eigv;
      eig = diag(3);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);
      }

    void eigenvaluetest()
      {
      mat matrix(4, 4), orig(4, 4);
      matrix(0, 0) = 4;
      matrix(0, 1) = 1;
      matrix(0, 2) = -2;
      matrix(0, 3) = 2;
      matrix(1, 0) = 1;
      matrix(1, 1) = 2;
      matrix(1, 2) = 0;
      matrix(1, 3) = 1;
      matrix(2, 0) = -2;
      matrix(2, 1) = 0;
      matrix(2, 2) = 3;
      matrix(2, 3) = -2;
      matrix(3, 0) = 2;
      matrix(3, 1) = 1;
      matrix(3, 2) = -2;
      matrix(3, 3) = -1;
      orig = matrix;
      mat eigs, eigv(4), result(4);
      eig_symm(eigs, matrix);
      eigv(0) = matrix(0, 0);
      eigv(1) = matrix(1, 0);
      eigv(2) = matrix(2, 0);
      eigv(3) = matrix(3, 0);
      result = orig * eigv;
      double eig = eigs(0);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);
      eigv(0) = matrix(0, 1);
      eigv(1) = matrix(1, 1);
      eigv(2) = matrix(2, 1);
      eigv(3) = matrix(3, 1);
      result = orig * eigv;
      eig = eigs(1);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);
      eigv(0) = matrix(0, 2);
      eigv(1) = matrix(1, 2);
      eigv(2) = matrix(2, 2);
      eigv(3) = matrix(3, 2);
      result = orig * eigv;
      eig = eigs(2);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);
      eigv(0) = matrix(0, 3);
      eigv(1) = matrix(1, 3);
      eigv(2) = matrix(2, 3);
      eigv(3) = matrix(3, 3);
      result = orig * eigv;
      eig = eigs(3);
      TEST_EQ_CLOSE(result(0), eig*eigv(0), 1e-8);
      TEST_EQ_CLOSE(result(1), eig*eigv(1), 1e-8);
      TEST_EQ_CLOSE(result(2), eig*eigv(2), 1e-8);
      TEST_EQ_CLOSE(result(3), eig*eigv(3), 1e-8);
      }


    void balanctest()
      {
      mat a(5, 5);
      a << 1, 2, 3, -7, 12, 2, 4, 7, 3, -1, 3, 7, 10, 8, 4, -7, 3, 8, -0.75, -9, 12, -1, 4, -9, 10;
      mat b = a;
      balanc(a);
      TEST_ASSERT(a == b);
      }

    void elmhestest()
      {
      mat a(5, 5);
      a << 1, 2, 3, -7, 12, 2, 4, 7, 3, -1, 3, 7, 10, 8, 4, -7, 3, 8, -0.75, -9, 12, -1, 4, -9, 10;
      elmhes(a);
      double tol = 1e-6;
      TEST_EQ_CLOSE(1.0, a(0, 0), tol);
      TEST_EQ_CLOSE(17.166667, a(0, 1), tol);
      TEST_EQ_CLOSE(-9.738494, a(0, 2), tol);
      TEST_EQ_CLOSE(2.674367, a(0, 3), tol);
      TEST_EQ_CLOSE(3.0, a(0, 4), tol);
      TEST_EQ_CLOSE(12.0, a(1, 0), tol);
      TEST_EQ_CLOSE(16.083333, a(1, 1), tol);
      TEST_EQ_CLOSE(-9.322176, a(1, 2), tol);
      TEST_EQ_CLOSE(-0.100844, a(1, 3), tol);
      TEST_EQ_CLOSE(4.0, a(1, 4), tol);
      TEST_EQ_CLOSE(0.25, a(2, 0), tol);
      TEST_EQ_CLOSE(3.319444, a(2, 1), tol);
      TEST_EQ_CLOSE(-11.372036, a(2, 2), tol);
      TEST_EQ_CLOSE(4.739487, a(2, 3), tol);
      TEST_EQ_CLOSE(10.333333, a(2, 4), tol);
      TEST_EQ_CLOSE(-0.583333, a(3, 0), tol);
      TEST_EQ_CLOSE(-0.307531, a(3, 1), tol);
      TEST_EQ_CLOSE(-11.556061, a(3, 2), tol);
      TEST_EQ_CLOSE(9.893547, a(3, 3), tol);
      TEST_EQ_CLOSE(15.715481, a(3, 4), tol);
      TEST_EQ_CLOSE(0.166667, a(4, 0), tol);
      TEST_EQ_CLOSE(-0.907950, a(4, 1), tol);
      TEST_EQ_CLOSE(0.224789, a(4, 2), tol);
      TEST_EQ_CLOSE(8.506681, a(4, 3), tol);
      TEST_EQ_CLOSE(8.645156, a(4, 4), tol);
      }

    void hqrtest()
      {
      mat a(5, 5), wr, wi;
      a << 1, 2, 3, -7, 12, 2, 4, 7, 3, -1, 3, 7, 10, 8, 4, -7, 3, 8, -0.75, -9, 12, -1, 4, -9, 10;
      elmhes(a);
      hqr(a, wr, wi);
      double tol = 1e-6;
      TEST_EQ_CLOSE(-10.486545, wr(0), tol);
      TEST_EQ_CLOSE(-7.774580, wr(1), tol);
      TEST_EQ_CLOSE(23.755955, wr(2), tol);
      TEST_EQ_CLOSE(18.291821, wr(3), tol);
      TEST_EQ_CLOSE(0.463350, wr(4), tol);
      TEST_EQ(0.0, wi(0));
      TEST_EQ(0.0, wi(1));
      TEST_EQ(0.0, wi(2));
      TEST_EQ(0.0, wi(3));
      TEST_EQ(0.0, wi(4));
      }

    void eigtest()
      {
      mat a(5, 5), wr, wi;
      a << 1, 2, 3, -7, 12, 2, 4, 7, 3, -1, 3, 7, 10, 8, 4, -7, 3, 8, -0.75, -9, 12, -1, 4, -9, 10;
      eig(a, wr, wi);
      double tol = 1e-6;
      TEST_EQ_CLOSE(-10.486545, wr(0), tol);
      TEST_EQ_CLOSE(-7.774580, wr(1), tol);
      TEST_EQ_CLOSE(23.755955, wr(2), tol);
      TEST_EQ_CLOSE(18.291821, wr(3), tol);
      TEST_EQ_CLOSE(0.463350, wr(4), tol);
      TEST_EQ(0.0, wi(0));
      TEST_EQ(0.0, wi(1));
      TEST_EQ(0.0, wi(2));
      TEST_EQ(0.0, wi(3));
      TEST_EQ(0.0, wi(4));
      }

    void aliasing()
      {
      mat m1(2, 3);
      mat m2(3, 2);

      m1(0, 0) = 2.0;
      m1(1, 0) = 6.0;
      m1(0, 1) = 3.0;
      m1(1, 1) = 2.0;
      m1(0, 2) = 7.0;
      m1(1, 2) = 1.0;

      m2(0, 0) = 1.0;
      m2(1, 0) = -1.0;
      m2(2, 0) = 2.0;
      m2(0, 1) = -3.0;
      m2(1, 1) = 4.0;
      m2(2, 1) = 9.0;

      m1 = m1 * m2;
      TEST_EQ(2, m1.rows());
      TEST_EQ(2, m1.cols());
      TEST_EQ(13.0, m1(0, 0));
      TEST_EQ(6.0, m1(1, 0));
      TEST_EQ(69.0, m1(0, 1));
      TEST_EQ(-1.0, m1(1, 1));

      m2 = transpose(m2);
      TEST_EQ(2, m2.rows());
      TEST_EQ(3, m2.cols());
      TEST_EQ(1.0, m2(0, 0));
      TEST_EQ(-1.0, m2(0, 1));
      TEST_EQ(2.0, m2(0, 2));
      TEST_EQ(-3.0, m2(1, 0));
      TEST_EQ(4.0, m2(1, 1));
      TEST_EQ(9.0, m2(1, 2));

      mat a1(2, 2);
      mat a2(2, 2);
      a1 << 1, 2, 3, 4;
      a2 << 1, 2, -1, -1;
      a1 += a1 * a2;
      TEST_EQ(2, a1.rows());
      TEST_EQ(2, a1.cols());
      TEST_EQ(0.0, a1(0, 0));
      TEST_EQ(2.0, a1(0, 1));
      TEST_EQ(2.0, a1(1, 0));
      TEST_EQ(6.0, a1(1, 1));

      mat b1(2, 2);
      mat b2(2, 2);
      b1 << 1, 2, 3, 4;
      b2 << 1, 2, -1, -1;
      b1 -= b1 * b2;
      TEST_EQ(2, b1.rows());
      TEST_EQ(2, b1.cols());
      TEST_EQ(2.0, b1(0, 0));
      TEST_EQ(2.0, b1(0, 1));
      TEST_EQ(4.0, b1(1, 0));
      TEST_EQ(2.0, b1(1, 1));
      }

    void noaliasing()
      {
      mat m1(2, 2);
      mat m2(2, 2);
      m1 << 1, 2, 3, 4;
      m2 << 1, 2, -1, -1;

      // wrong results, because user requested noalias()
      m1.noalias() = m1 * m2;
      TEST_EQ(2, m1.rows());
      TEST_EQ(2, m1.cols());
      TEST_EQ(-1.0, m1(0, 0));
      TEST_EQ(-1.0, m1(1, 0));
      TEST_EQ(-4.0, m1(0, 1));
      TEST_EQ(-6.0, m1(1, 1));

      // wrong results, because user requested noalias()
      m2.noalias() = transpose(m2);
      TEST_EQ(2, m2.rows());
      TEST_EQ(2, m2.cols());
      TEST_EQ(1.0, m2(0, 0));
      TEST_EQ(-1.0, m2(0, 1));
      TEST_EQ(-1.0, m2(0, 2));
      TEST_EQ(-1.0, m2(1, 0));

      // wrong results, because user requested noalias()
      mat a1(2, 2);
      mat a2(2, 2);
      a1 << 1, 2, 3, 4;
      a2 << 1, 2, -1, -1;
      a1.noalias() += a1 * a2;
      TEST_EQ(2, a1.rows());
      TEST_EQ(2, a1.cols());
      TEST_EQ(0.0, a1(0, 0));
      TEST_EQ(0.0, a1(0, 1));
      TEST_EQ(2.0, a1(1, 0));
      TEST_EQ(4.0, a1(1, 1));

      // wrong results, because user requested noalias()
      mat b1(2, 2);
      mat b2(2, 2);
      b1 << 1, 2, 3, 4;
      b2 << 1, 2, -1, -1;
      b1.noalias() -= b1 * b2;
      TEST_EQ(2, b1.rows());
      TEST_EQ(2, b1.cols());
      TEST_EQ(2.0, b1(0, 0));
      TEST_EQ(0.0, b1(0, 1));
      TEST_EQ(4.0, b1(1, 0));
      TEST_EQ(0.0, b1(1, 1));
      }

    void datatest()
      {
      mat4 m(2, 2);
      m << 1, 2, 3, 4;
      double* p = m.data();
      TEST_EQ(1.0, *p++);
      TEST_EQ(2.0, *p++);
      TEST_EQ(3.0, *p++);
      TEST_EQ(4.0, *p++);
      }

    void tracetest()
      {
      matrix<double> m(5, 5);
      for (int i = 0; i < 5; ++i)
        {
        for (int j = 0; j < 5; ++j)
          {
          m(i, j) = (double)(i * 5 + j);
          }
        }
      m(0, 0) = 100.0;
      double tr = trace(m);
      TEST_EQ(100. + (1.*5. + 1.) + (2.*5. + 2.) + (3.*5. + 3.) + (4.*5. + 4.), tr);
      tr = trace(m + m);
      TEST_EQ((100. + (1.*5. + 1.) + (2.*5. + 2.) + (3.*5. + 3.) + (4.*5. + 4.))*2.0, tr);
      }

    void vconcattest()
      {
      matrix<double> m_up(3, 3);
      m_up << 1, 2, 3, 4, 5, 6, 7, 8, 9;
      matrix<double> m_mid(4, 3);
      m_mid << 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
      matrix<double> m_down(3, 3);
      m_down << 22, 23, 24, 25, 26, 27, 28, 29, 30;
      matrix<double> m_concat(10, 3);
      m_concat << m_up, m_mid, m_down;
      auto it = m_concat.begin();
      for (int i = 1; i < 31; ++i)
        {
        TEST_EQ(*it++, i);
        }
      TEST_ASSERT(it == m_concat.end());
      }

    void vconcattest2()
      {
      matrix<double> m_1(9, 3);
      int i = 1;
      for (auto& v : m_1)
        v = i++;

      matrix<double> m_2(9, 3);

      m_2 << block(m_1, 0, 0, 3, 3), block(m_1, 3, 0, 3, 3), block(m_1, 6, 0, 3, 3);

      TEST_ASSERT(m_1 == m_2);
      }

    void hconcattest()
      {
      matrix<double> m_left(3, 3);
      m_left << 1, 2, 3, 4, 5, 6, 7, 8, 9;
      matrix<double> m_mid(3, 4);
      m_mid << 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
      matrix<double> m_right(3, 3);
      m_right << 22, 23, 24, 25, 26, 27, 28, 29, 30;
      matrix<double> m_concat(3, 10);
      m_concat << m_left, m_mid, m_right;
      auto it = m_concat.begin();
      for (int r = 0; r < 3; ++r)
        {
        for (int c = 0; c < 3; ++c)
          {
          TEST_EQ(m_concat(r, c), m_left(r, c));
          TEST_EQ(m_concat(r, c + 3), m_mid(r, c));
          TEST_EQ(m_concat(r, c + 7), m_right(r, c));
          }
        }
      for (int r = 0; r < 3; ++r)
        {
        TEST_EQ(m_concat(r, 6), m_mid(r, 3));
        }
      }

    void hconcattest2()
      {
      matrix<double> m_1(3, 9);
      int i = 1;
      for (auto& v : m_1)
        v = i++;

      matrix<double> m_2(3, 9);

      m_2 << block(m_1, 0, 0, 3, 3), block(m_1, 0, 3, 3, 3), block(m_1, 0, 6, 3, 3);

      TEST_ASSERT(m_1 == m_2);
      }

    struct sparse_vector_fixture
      {
      sparse_vector_fixture() : vector_1(5), vector_2(5), vector_3(10), vector_4(5)
        {
        size_t i;
        for (i = 0; i < 5; ++i)
          {
          vector_1.put(i) = static_cast<double>(i);
          vector_2.put(i) = static_cast<double>(i)*2.0;
          vector_4.put(i) = static_cast<int>(i) * 4;
          }
        for (i = 0; i < 10; ++i)
          {
          vector_3.put(i) = static_cast<double>(i)*3.0;
          }
        }

      sparse_vector<double> vector_1;
      sparse_vector<double> vector_2;
      sparse_vector<double> vector_3;
      sparse_vector<int> vector_4;

      };

    void sparse_vector_construction()
      {
      sparse_vector<double> v(5);
      size_t i;
      for (i = 0; i < 5; ++i)
        v.put(i) = static_cast<double>(i);
      for (i = 0; i < 5; ++i)
        TEST_EQ(v.get(i), i);
      sparse_vector<double> v2(0);
      TEST_ASSERT(v2.begin() == v2.end());
      }

    void sparse_vector_sparsity_optimization()
      {
      sparse_vector<double> v(5);
      size_t i;
      for (i = 0; i < 5; ++i)
        v.put(i) = 0.0;
      v.put(0) = 1.0;
      TEST_EQ(v.entries_stored(), static_cast<size_t>(5));
      sparse_vector<double> v2 = v;
      TEST_EQ(v2.entries_stored(), static_cast<size_t>(1));
      TEST_EQ(v2.get(0), 1.0);
      }

    struct sparse_vector_copy_constructor : public sparse_vector_fixture
      {
      void test()
        {
        int i;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), i);
        sparse_vector<double> newcontainer(vector_1);
        for (i = 0; i < 5; ++i)
          TEST_EQ(newcontainer.get(i), i);
        }
      };

    struct sparse_vector_swap : public sparse_vector_fixture
      {
      void test()
        {
        vector_1.swap(vector_3);
        int i;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_3.get(i), i);
        for (i = 0; i < 10; ++i)
          TEST_EQ(vector_1.get(i), i * 3);
        }
      };

    struct sparse_vector_assignment : public sparse_vector_fixture
      {
      void test()
        {
        vector_3 = vector_1;
        size_t i;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_3.get(i), i);
        TEST_EQ(vector_3.size(), static_cast<size_t>(5));
        }
      };

    struct sparse_vector_mul_div : public sparse_vector_fixture
      {
      void test()
        {
        vector_1 *= 2.0;
        size_t i;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), 2.0*i);
        vector_1 /= 2.0;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), i);
        }
      };

    struct sparse_vector_add : public sparse_vector_fixture
      {
      void test()
        {
        vector_1 += vector_2;
        size_t i;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), 3.0*i);
        vector_1 -= vector_2;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), i);
        vector_1 += vector_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), 2.0*i);
        vector_1 -= vector_1;
        vector_1 -= vector_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(vector_1.get(i), 0.0);
        }
      };

    struct sparse_vector_equality : public sparse_vector_fixture
      {
      void test()
        {
        TEST_ASSERT(vector_1 == vector_1);
        TEST_ASSERT(vector_1 != vector_2);
        }
      };

    struct fixture_sparse_matrix
      {
      fixture_sparse_matrix() : matrix_1(3, 3), matrix_2(3, 3), matrix_3(4, 1), matrix_4(1, 2)
        {
        size_t i, j;
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            matrix_1.put(i, j) = static_cast<double>(i * 3 + j)*1.0;
            matrix_2.put(i, j) = static_cast<double>(i * 3 + j)*2.0;
            }
        for (i = 0; i < 4; ++i)
          {
          matrix_3.put(i, 0) = static_cast<double>(i)*3.0;
          }
        }


      sparse_matrix<double> matrix_1;
      sparse_matrix<double> matrix_2;
      sparse_matrix<double> matrix_3;
      sparse_matrix<double> matrix_4;
      };

    struct sparse_matrix_getput_test : public fixture_sparse_matrix
      {
      void test()
        {
        int i, j;
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_1.get(i, j), static_cast<double>(i * 3 + j));
            }
        matrix_1.put(1, 1) = 100.0;
        TEST_EQ(matrix_1.get(1, 1), 100.0);
        matrix_1.put(1, 2) = 101.0;
        TEST_EQ(matrix_1.get(1, 2), 101.0);
        }
      };

    struct sparse_matrix_size_test : public fixture_sparse_matrix
      {
      void test()
        {
        TEST_EQ(matrix_1.rows(), static_cast<size_t>(3));
        TEST_EQ(matrix_1.cols(), static_cast<size_t>(3));

        TEST_EQ(matrix_3.rows(), static_cast<size_t>(4));
        TEST_EQ(matrix_3.cols(), static_cast<size_t>(1));

        TEST_EQ(matrix_4.rows(), static_cast<size_t>(1));
        TEST_EQ(matrix_4.cols(), static_cast<size_t>(2));
        }
      };

    struct sparse_matrix_swap_test : public fixture_sparse_matrix
      {
      void test()
        {
        int i, j;
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_1.get(i, j), static_cast<double>(i * 3 + j));
            TEST_EQ(matrix_2.get(i, j), static_cast<double>(i * 3 + j) * 2);
            }
        matrix_1.swap(matrix_2);
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_2.get(i, j), static_cast<double>(i * 3 + j));
            TEST_EQ(matrix_1.get(i, j), static_cast<double>(i * 3 + j) * 2);
            }
        matrix_1.swap(matrix_4);
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_2.get(i, j), static_cast<double>(i * 3 + j));
            TEST_EQ(matrix_4.get(i, j), static_cast<double>(i * 3 + j) * 2);
            }
        TEST_EQ(matrix_1.rows(), static_cast<size_t>(1));
        TEST_EQ(matrix_1.cols(), static_cast<size_t>(2));
        }
      };

    struct sparse_matrix_assignment_test : public fixture_sparse_matrix
      {
      void test()
        {
        int i, j;
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_1.get(i, j), static_cast<double>(i * 3 + j));
            TEST_EQ(matrix_2.get(i, j), static_cast<double>(i * 3 + j) * 2);
            }
        matrix_1 = matrix_2;
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_2.get(i, j), static_cast<double>(i * 3 + j) * 2);
            TEST_EQ(matrix_1.get(i, j), static_cast<double>(i * 3 + j) * 2);
            }
        matrix_1 = matrix_4;
        for (i = 0; i < 3; ++i)
          for (j = 0; j < 3; ++j)
            {
            TEST_EQ(matrix_2.get(i, j), static_cast<double>(i * 3 + j) * 2);
            }
        TEST_EQ(matrix_1.rows(), static_cast<size_t>(1));
        TEST_EQ(matrix_1.cols(), static_cast<size_t>(2));

        TEST_EQ(matrix_4.rows(), static_cast<size_t>(1));
        TEST_EQ(matrix_4.cols(), static_cast<size_t>(2));
        }
      };

    struct sparse_matrix_resize_test : public fixture_sparse_matrix
      {
      void test()
        {
        matrix_1.resize(100, 500);
        TEST_EQ(matrix_1.rows(), static_cast<size_t>(100));
        TEST_EQ(matrix_1.cols(), static_cast<size_t>(500));
        }
      };

    struct sparse_matrix_iterator_test : public fixture_sparse_matrix
      {
      void test()
        {
        auto iter = matrix_1.begin();
        auto end = matrix_1.end();
        int i = 0;
        for (; iter != end; ++iter)
          {
          TEST_EQ((*iter), static_cast<double>(i++));
          }
        auto citer = matrix_1.cbegin();
        auto cend = matrix_1.cend();
        i = 0;
        for (; citer != cend; ++citer)
          {
          TEST_EQ((*citer), static_cast<double>(i++));
          }
        iter = matrix_1.begin();
        (*iter) = 5.0;
        TEST_EQ(matrix_1(0, 0), 5.0);
        }
      };

    struct fixture_sparse_matrix_operations
      {
      fixture_sparse_matrix_operations() : mat_1(5, 5), mat_2(5, 5), smat_1(5, 5), smat_2(5, 5), mat_1col(5, 5)
        {
        size_t i, j;
        for (i = 0; i < 5; ++i)
          for (j = 0; j < 5; ++j)
            {
            mat_1(i, j) = (double)(i * 5 + j);
            mat_2(i, j) = (double)((i * 5 + j)*2.0);
            mat_1col(i, j) = (double)(i + 5 * j);
            }
        for (i = 0; i < 5; ++i)
          {
          smat_1.put(i, i) = (double)(i);
          smat_2.put(i, i) = (double)(i*2.0);

          }
        }

      matrix<double> mat_1;
      matrix<double> mat_2;
      matrix<double> mat_1col;
      sparse_matrix<double> smat_1;
      sparse_matrix<double> smat_2;
      };

    struct sparse_matrix_add_test : public fixture_sparse_matrix_operations
      {
      void test()
        {
        int i;
        sparse_matrix<double> sm = smat_1 + smat_2;

        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*3.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = smat_1 + (smat_1 + smat_2);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*4.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = (smat_1 + smat_2) + smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*4.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = (smat_1 + smat_2) + (smat_2 + smat_1);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*6.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = smat_1 + smat_2 + smat_2 + smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*6.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        }
      };

    struct sparse_matrix_subtract_test : public fixture_sparse_matrix_operations
      {
      void test()
        {
        int i;
        sparse_matrix<double> sm = smat_1 - smat_2;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), -static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = smat_1 - (smat_1 + smat_2);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), -static_cast<double>(i*2.0));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = (smat_1 + smat_2) - smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*2.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = (smat_1 + smat_2) - (smat_2 - smat_1);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), i*2.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = smat_1 - smat_2 - smat_2 - smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), -static_cast<double>(i*4.0));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        }
      };

    struct sparse_matrix_negate_test : public fixture_sparse_matrix_operations
      {
      void test()
        {
        int i;
        sparse_matrix<double> sm = -smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), -static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = -(smat_1 + smat_2);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), -static_cast<double>(3.0*i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));
        }
      };

    struct sparse_matrix_scalar_mul_test : public fixture_sparse_matrix_operations
      {
      void test()
        {
        int i;
        sparse_matrix<double> sm = smat_1 * 2.0;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 2.0*static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = 2.0 * smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 2.0*static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = (smat_1 + smat_2)*2.0;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 6.0*static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = 2.0*(smat_1 + smat_2);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 6.0*static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));
        }
      };

    struct sparse_matrix_scalar_div_test : public fixture_sparse_matrix_operations
      {
      void test()
        {
        int i;
        sparse_matrix<double> sm = 2.0 / smat_1;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 2.0 / static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(5));

        sm = smat_1 / 2.0;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), static_cast<double>(i) / 2.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = (smat_1 + smat_2) / 2.0;
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 3.0*static_cast<double>(i) / 2.0);
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));

        sm = 2.0 / (smat_1 + smat_2);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 2.0 / (3.0*static_cast<double>(i)));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(5));
        }
      };

    struct sparse_matrix_transpose_test : public fixture_sparse_matrix_operations
      {
      void test()
        {
        int i;
        sparse_matrix<double> sm = transpose(smat_1);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(5));

        sm = transpose(smat_1 + smat_2);
        for (i = 0; i < 5; ++i)
          TEST_EQ(sm.get(i, i), 3.0*static_cast<double>(i));
        TEST_EQ(sm.entries_stored(), static_cast<size_t>(4));
        }
      };

    void sparse_matrix_multiply()
      {
      sparse_matrix<double> sm1(2, 3);
      sm1.put(0, 0) = 1.0;
      sm1.put(0, 2) = 2.0;
      sm1.put(1, 0) = -1.0;
      sm1.put(1, 1) = 3.0;
      sm1.put(1, 2) = 1.0;
      TEST_EQ(sm1.entries_stored(), static_cast<size_t>(5));

      sparse_matrix<double> sm2(3, 2);
      sm2.put(0, 0) = 3.0;
      sm2.put(0, 1) = 1.0;
      sm2.put(1, 0) = 2.0;
      sm2.put(1, 1) = 1.0;
      sm2.put(2, 0) = 1.0;
      TEST_EQ(sm2.entries_stored(), static_cast<size_t>(5));

      sparse_matrix<double> sm = sm1 * sm2;
      TEST_EQ(sm.get(0, 0), 5.0);
      TEST_EQ(sm.get(0, 1), 1.0);
      TEST_EQ(sm.get(1, 0), 4.0);
      TEST_EQ(sm.get(1, 1), 2.0);
      TEST_EQ(sm.rows(), static_cast<uint64_t>(2));
      TEST_EQ(sm.cols(), static_cast<uint64_t>(2));
      TEST_EQ(sm.entries_stored(), static_cast<uint64_t>(4));

      sm = (sm1 + sm1)*sm2;
      TEST_EQ(sm.get(0, 0), 2.0*5.0);
      TEST_EQ(sm.get(0, 1), 2.0*1.0);
      TEST_EQ(sm.get(1, 0), 2.0*4.0);
      TEST_EQ(sm.get(1, 1), 2.0*2.0);
      TEST_EQ(sm.rows(), static_cast<size_t>(2));
      TEST_EQ(sm.cols(), static_cast<size_t>(2));

      sm = sm1 * (sm2 + sm2);
      TEST_EQ(sm.get(0, 0), 2.0*5.0);
      TEST_EQ(sm.get(0, 1), 2.0*1.0);
      TEST_EQ(sm.get(1, 0), 2.0*4.0);
      TEST_EQ(sm.get(1, 1), 2.0*2.0);
      TEST_EQ(sm.rows(), static_cast<size_t>(2));
      TEST_EQ(sm.cols(), static_cast<size_t>(2));

      sm = (2.0*sm1)*(sm2 + sm2);
      TEST_EQ(sm.get(0, 0), 4.0*5.0);
      TEST_EQ(sm.get(0, 1), 4.0*1.0);
      TEST_EQ(sm.get(1, 0), 4.0*4.0);
      TEST_EQ(sm.get(1, 1), 4.0*2.0);
      TEST_EQ(sm.rows(), static_cast<size_t>(2));
      TEST_EQ(sm.cols(), static_cast<size_t>(2));
      }

    void sparse_matrix_vector_multiply()
      {
      sparse_matrix<double> sm1(2, 3);
      sm1.put(0, 0) = 1.0;
      sm1.put(0, 2) = 2.0;
      sm1.put(1, 0) = -1.0;
      sm1.put(1, 1) = 3.0;
      sm1.put(1, 2) = 1.0;

      mat v1(3);
      v1(0) = 3.0;
      v1(1) = 2.0;
      v1(2) = 1.0;

      mat x = sm1 * v1;
      TEST_EQ(x(0), 5.0);
      TEST_EQ(x(1), 4.0);
      TEST_EQ(x.rows(), static_cast<size_t>(2));
      TEST_EQ(x.cols(), static_cast<size_t>(1));

      x = (sm1 + sm1) * v1;
      TEST_EQ(x(0), 10.0);
      TEST_EQ(x(1), 8.0);
      TEST_EQ(x.rows(), static_cast<size_t>(2));
      TEST_EQ(x.cols(), static_cast<size_t>(1));

      x = sm1 * (v1 + v1);
      TEST_EQ(x(0), 10.0);
      TEST_EQ(x(1), 8.0);
      TEST_EQ(x.rows(), static_cast<size_t>(2));
      TEST_EQ(x.cols(), static_cast<size_t>(1));

      x = (sm1 + sm1) * (v1 + v1);
      TEST_EQ(x(0), 20.0);
      TEST_EQ(x(1), 16.0);
      TEST_EQ(x.rows(), static_cast<size_t>(2));
      TEST_EQ(x.cols(), static_cast<size_t>(1));
      }

    void sparse_norm_tests()
      {
      sparse_matrix<double> v(5, 1);
      v.put(3, 1) = 2.0;
      v.put(4, 1) = 3.0;
      TEST_EQ(std::sqrt(4 + 9), norm(v));
      sparse_matrix<double> m(2, 2);
      m.put(0, 1) = 5.0;
      m.put(1, 1) = 7.0;
      TEST_EQ(std::sqrt(25 + 49), norm(m));
      sparse_matrix<double> m2(2, 2);
      m2.put(0, 0) = -2;
      TEST_EQ(std::sqrt(4 + 25 + 49), norm(m + m2));
      }

    void sparse_diagonal_tests()
      {
      sparse_matrix<double> m(4, 5);
      m.put(0, 0) = 5.0;
      m.put(1, 1) = 7.0;
      m.put(3, 1) = 17.0;
      m.put(0, 2) = 9.0;
      mat d = diagonal(m);
      TEST_EQ(d.rows(), 4);
      TEST_EQ(d.cols(), 1);
      TEST_EQ(5.0, d(0));
      TEST_EQ(7.0, d(1));
      TEST_EQ(0.0, d(2));
      TEST_EQ(0.0, d(3));
      d = diagonal(m + m);
      TEST_EQ(10.0, d(0));
      TEST_EQ(14.0, d(1));
      TEST_EQ(0.0, d(2));
      TEST_EQ(0.0, d(3));
      }

    void sparse_block_matrix()
      {
      sparse_matrix<double> m(4, 4);
      for (int i = 0; i < 16; ++i)
        m.put(i / 4, i % 4) = i + 1;

      sparse_matrix<double> b = block(m, 1, 1, 2, 2);
      TEST_EQ(6.f, b(0, 0));
      TEST_EQ(7.f, b(0, 1));
      TEST_EQ(10.f, b(1, 0));
      TEST_EQ(11.f, b(1, 1));

      b = block(m, 0, 0, 3, 3);
      TEST_EQ(1.f, b(0, 0));
      TEST_EQ(2.f, b(0, 1));
      TEST_EQ(3.f, b(0, 2));
      TEST_EQ(5, b(1, 0));
      TEST_EQ(6, b(1, 1));
      TEST_EQ(7, b(1, 2));
      TEST_EQ(9, b(2, 0));
      TEST_EQ(10.f, b(2, 1));
      TEST_EQ(11.f, b(2, 2));

      mat d = diagonal(block(m, 0, 0, 3, 3));
      TEST_EQ(1.f, d(0));
      TEST_EQ(6.f, d(1));
      TEST_EQ(11.f, d(2));

      b = block(block(m, 0, 0, 3, 3), 1, 1, 2, 2);
      TEST_EQ(6.f, b(0, 0));
      TEST_EQ(7.f, b(0, 1));
      TEST_EQ(10.f, b(1, 0));
      TEST_EQ(11.f, b(1, 1));
      }

    void sparse_trace_test()
      {
      sparse_matrix<double> m(5, 5);
      for (int i = 0; i < 5; ++i)
        {
        for (int j = 0; j < 5; ++j)
          {
          m.put(i, j) = (double)(i * 5 + j);
          }
        }
      m.put(0, 0) = 100.0;
      double tr = trace(m);
      TEST_EQ(100. + (1.*5. + 1.) + (2.*5. + 2.) + (3.*5. + 3.) + (4.*5. + 4.), tr);
      tr = trace(m + m);
      TEST_EQ((100. + (1.*5. + 1.) + (2.*5. + 2.) + (3.*5. + 3.) + (4.*5. + 4.))*2.0, tr);
      }

    void sparse_matrix_init()
      {
      sparse_matrix<double> m = zeros(3, 4);
      TEST_EQ(3, m.rows());
      TEST_EQ(4, m.cols());
      for (int r = 0; r < m.rows(); ++r)
        {
        for (int c = 0; c < m.cols(); ++c)
          TEST_EQ(0.f, m(r, c));
        }
      TEST_EQ(0, m.entries_stored());
      
      m = ones(5, 7);
      TEST_EQ(5, m.rows());
      TEST_EQ(7, m.cols());
      for (int r = 0; r < m.rows(); ++r)
        {
        for (int c = 0; c < m.cols(); ++c)
          TEST_EQ(1.f, m(r, c));
        }
      
      m = identity(2, 3);
      TEST_EQ(2, m.rows());
      TEST_EQ(3, m.cols());
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, m(i, j));      
      TEST_EQ(2, m.entries_stored());
      }

    void assign_sparse_to_dense()
      {
      sparse_matrix<double> m = identity(3, 4);
      matrix<double> d = m;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, d(i, j));

      m = d;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
          TEST_EQ(i == j ? 1.0 : 0.0, m(i, j));
      TEST_EQ(3, m.entries_stored());
      }

    void conjugate_gradient_tests()
      {
      smat A(10, 10);
      for (int i = 1; i <= 10; ++i)
        A.put(i - 1, i - 1) = (double)i;
      mat b(10);
      for (int i = 1; i <= 10; ++i)
        b(i - 1) = (double)i;

      mat x;
      mat x0 = zeros(10, 1);
      double residu;
      uint64_t iter;
      conjugate_gradient(x, residu, iter, A, b, x0, 0.00001);

      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0, x(i), 0.00001);
      TEST_EQ(9, iter);

      b = ones(10, 1);
      conjugate_gradient(x, residu, iter, A, b, x0, 0.00001);
      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0 / ((double)(i + 1)), x(i), 0.00001);
      TEST_EQ(9, iter);
      }

    void preconditioned_conjugate_gradient_tests()
      {
      smat A(10, 10);
      smat P(10, 10);
      for (int i = 1; i <= 10; ++i)
        {
        A.put(i - 1, i - 1) = (double)i;
        P.put(i - 1, i - 1) = 1.0/(double)i;
        }

      mat b(10);
      for (int i = 1; i <= 10; ++i)
        b(i - 1) = (double)i;

      mat x;
      mat x0 = zeros(10, 1);
      double residu;
      uint64_t iter;
      preconditioned_conjugate_gradient(x, residu, iter, A, P, b, x0, 0.00001);

      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0, x(i), 0.00001);
      TEST_EQ(1, iter);

      b = ones(10, 1);
      preconditioned_conjugate_gradient(x, residu, iter, A, P, b, x0, 0.00001);
      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0 / ((double)(i + 1)), x(i), 0.00001);
      TEST_EQ(1, iter);
      }

    void bicgstab_tests()
      {
      smat A(10, 10);
      for (int i = 1; i <= 10; ++i)
        A.put(i - 1, i - 1) = (double)i;
      mat b(10);
      for (int i = 1; i <= 10; ++i)
        b(i - 1) = (double)i;

      mat x;
      mat x0 = zeros(10, 1);
      double residu;
      uint64_t iter;
      bicgstab(x, residu, iter, A, b, x0, 0.000001);

      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0, x(i), 0.00001);
      TEST_EQ(10, iter);
      
      b = ones(10, 1);
      bicgstab(x, residu, iter, A, b, x0, 0.000001);
      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0 / ((double)(i + 1)), x(i), 0.00001);
      TEST_EQ(9, iter);      
      }

    void bipcgstab_tests()
      {
      smat A(10, 10);
      smat P(10, 10);
      for (int i = 1; i <= 10; ++i)
        {
        A.put(i - 1, i - 1) = (double)i;
        P.put(i - 1, i - 1) = 1.0 / (double)i;
        }
      mat b(10);
      for (int i = 1; i <= 10; ++i)
        b(i - 1) = (double)i;

      mat x;
      mat x0 = zeros(10, 1);
      double residu;
      uint64_t iter;
      bipcgstab(x, residu, iter, A, P, b, x0, 0.000001);

      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0, x(i), 0.00001);
      TEST_EQ(1, iter);

      b = ones(10, 1);
      bipcgstab(x, residu, iter, A, P, b, x0, 0.000001);
      for (int i = 0; i < 10; ++i)
        TEST_EQ_CLOSE(1.0 / ((double)(i + 1)), x(i), 0.00001);
      TEST_EQ(1, iter);
      }      
    }
  }

namespace
  {
  void negating_outside_namespace()
    {
    jtk::mat m(3, 3);
    m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    jtk::mat m2 = -jtk::transpose(m);
    TEST_EQ(-1.0, m2(0, 0));
    TEST_EQ(-2.0, m2(1, 0));
    TEST_EQ(-3.0, m2(2, 0));
    TEST_EQ(-4.0, m2(0, 1));
    TEST_EQ(-5.0, m2(1, 1));
    TEST_EQ(-6.0, m2(2, 1));
    TEST_EQ(-7.0, m2(0, 2));
    TEST_EQ(-8.0, m2(1, 2));
    TEST_EQ(-9.0, m2(2, 2));
    }
  }


void run_all_mat_tests()
  {
  using namespace jtk;
  matrix_construction_vector().test();
  matrix_construction_array().test();
  copy_constructor().test();
  copy_constructor_array().test();
  copy_assignment_constructor().test();
  copy_assignment_array().test();
  vector_test();
  vector_test_array();
  svd_1();
  svd_2();
  svd_3();
  pseudo_inverse_1();
  pseudo_inverse_2();
  pseudo_inverse_3();
  svd_1_array();
  svd_2_array();
  svd_3_array();
  pseudo_inverse_1_array();
  pseudo_inverse_2_array();
  pseudo_inverse_3_array();
  adding().test();
  adding2().test();
  adding3().test();
  adding4().test();
  adding5().test();
  adding6().test();
  subtracting1().test();
  subtracting2().test();
  subtracting3().test();
  subtracting4().test();
  negating1().test();
  negating2().test();
  negating3().test();
  scalarmultiply().test();
  scalardivision();
  scalaraddition();
  scalarsubtraction();
  matrix_multiplication();
  adding_array().test();
  adding2_array().test();
  adding3_array().test();
  adding4_array().test();
  adding5_array().test();
  adding6_array().test();
  subtracting1_array().test();
  subtracting2_array().test();
  subtracting3_array().test();
  subtracting4_array().test();
  negating1_array().test();
  negating2_array().test();
  scalarmultiply_array().test();
  scalardivision_array();
  scalaraddition_array();
  scalarsubtraction_array();
  matrix_multiplication_array();
  matrixconversion().test();
  matrixtranspose().test();
  matrixtranspose_array().test();
  lsd_test();
  lsd_test_array();
  outputstream();
  inputstream();
  matrix_add_assignment().test();
  matrix_subtract_assignment().test();
  matrix_multiply_assignment();
  matrix_multiply_assignment2();
  matrix_add_asignment_array().test();
  matrix_subtract_assignment_array().test();
  matrix_multiply_assignment_array();
  matrix_multiply_assignment2_array();
  matrix_diagonal().test();
  block_matrix();
  matrix_comparison().test();
  matrix_capacity().test();
  matrix_init_float();
  matrix_init_double();
  lu_dcmp();
  lu_dcmp2();
  lu_bksb();
  lu_solve();
  lu_invert();
  lu_determinant();
  lu_dcmp_array();
  lu_dcmp2_array();
  lu_bksb_array();
  lu_solve_array();
  lu_invert_array();
  lu_determinant_array();
  cholesky1();
  back_substitution();
  forward_substitution();
  solve_cholesky();
  solve_cholesky_2();
  cholesky1_array();
  cholesky2_array();
  back_substitution_array();
  forward_substitution_array();
  solve_cholesky_array();
  solve_cholesky_2_array();
  very_large_size_expression_uncomputed();
  solve_qr();
  qr();
  qr_rectangular();
  qr_rectangular2();
  solve_qr_test_least_squares();
  solve_qr_test_least_squares_2();
  solve_qr_array();
  qr_array();
  qr_rectangular_array();
  qr_rectangular2_array();
  solve_qr_test_least_squares_array();
  solve_qr_test_least_squares_2_array();
  norm_tests();
  blocknorm();
  norm_tests_array();
  blocknorm_array();
  qrfac();
  qrfac2();
  qrsolvtest();
  qrsolvtest2();
  qrsolvtest3();
  qrfac_array();
  qrfac2_array();
  qrsolvtest_array();
  fdjac_test();
  fdjac_test_array();
  lm_test();
  lm_test_array();
  tred2test();
  tqlitest();
  eigenvaluetest();
  balanctest();
  elmhestest();
  hqrtest();
  eigtest();
  negating_outside_namespace();
  transpose_aliasing_problem();
  transpose_aliasing_problem_2();
  aliasing();
  noaliasing();
  datatest();
  tracetest();
  vconcattest();
  hconcattest();
  vconcattest2();
  hconcattest2();
  sparse_vector_construction();
  sparse_vector_sparsity_optimization();
  sparse_vector_copy_constructor().test();
  sparse_vector_swap().test();
  sparse_vector_assignment().test();
  sparse_vector_mul_div().test();
  sparse_vector_add().test();
  sparse_vector_equality().test();  
  sparse_matrix_getput_test().test();
  sparse_matrix_size_test().test();
  sparse_matrix_swap_test().test();
  sparse_matrix_assignment_test().test();
  sparse_matrix_resize_test().test();
  sparse_matrix_iterator_test().test();
  sparse_matrix_add_test().test();
  sparse_matrix_subtract_test().test();
  sparse_matrix_negate_test().test();
  sparse_matrix_scalar_mul_test().test();
  sparse_matrix_scalar_div_test().test();
  sparse_matrix_transpose_test().test();
  sparse_matrix_multiply();
  sparse_matrix_vector_multiply();
  sparse_norm_tests();
  sparse_diagonal_tests();
  sparse_block_matrix();
  sparse_trace_test();
  sparse_matrix_init();
  assign_sparse_to_dense();
  conjugate_gradient_tests();
  preconditioned_conjugate_gradient_tests();
  bicgstab_tests();
  bipcgstab_tests();  
  }