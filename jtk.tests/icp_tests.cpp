#include "icp_tests.h"
#include "test_assert.h"

#define JTK_ICP_IMPLEMENTATION
#include "../jtk/icp.h"
#include <random>

namespace jtk
  {
  namespace
    {
    template <class T>
    T get_random_number(T low = 0, T high = 1)
      {
      static std::default_random_engine generator = std::default_random_engine();
      std::uniform_real_distribution<T> dist(low, high);
      return dist(generator);
      }

    int get_random_integer(int low, int high)
      {
      static std::default_random_engine generator = std::default_random_engine();
      std::uniform_int_distribution<int> dist(low, high);
      return dist(generator);
      }

    template <class T>
    jtk::vec3<T> get_random_unit_vector_3()
      {
      jtk::vec3<T> v;
      for (int i = 0; i < 3; ++i)
        v[i] = get_random_number<T>(0, 1);

      auto l = jtk::length_sqr(v);

      if (l < 1e-10)
        return get_random_unit_vector_3<T>();

      return v / std::sqrt(l);
      }


    void test_icp_1()
      {
      std::vector<jtk::vec3<float>> pts_1, pts_2;

      for (int i = 0; i < 100; ++i)
        {
        pts_1.push_back(get_random_unit_vector_3<float>());
        pts_2.push_back(pts_1.back() + jtk::vec3<float>(1, 2, 3));
        }
      matf16 t = identity<float>(4,4);
      t(0, 3) = 0.9f;
      t(1, 3) = 2.3f;
      t(2, 3) = 2.7f;      
      float residual;
      auto transf = icp(residual, pts_2, pts_1, t);
      matf16 exact_solution = identity<float>(4, 4);
      exact_solution(0, 3) = 1.f;
      exact_solution(1, 3) = 2.f;
      exact_solution(2, 3) = 3.f;
      for (int r = 0; r < 4; ++r)
        {
        for (int c = 0; c < 4; ++c)
          {
          TEST_EQ_CLOSE(exact_solution(r, c), transf(r, c), 1e-3);
          }
        }
      TEST_ASSERT(residual < 0.001f);
      }
    }

  void test_icp_2()
    {
    std::vector<jtk::vec3<double>> pts_1, pts_2;

    for (int i = 0; i < 100; ++i)
      {
      pts_1.push_back(get_random_unit_vector_3<double>());
      pts_2.push_back(pts_1.back() + jtk::vec3<double>(1, 2, 3));
      }
    mat16 t = identity<double>(4, 4);
    t(0, 3) = 0.9;
    t(1, 3) = 2.3;
    t(2, 3) = 2.7;
    double residual;
    auto transf = icp(residual, pts_2, pts_1, t);
    mat16 exact_solution = identity<double>(4, 4);
    exact_solution(0, 3) = 1;
    exact_solution(1, 3) = 2;
    exact_solution(2, 3) = 3;
    for (int r = 0; r < 4; ++r)
      {
      for (int c = 0; c < 4; ++c)
        {
        TEST_EQ_CLOSE(exact_solution(r, c), transf(r, c), 1e-3);
        }
      }
    TEST_ASSERT(residual < 0.001);
    }  
  }

void run_all_icp_tests()
  {
  using namespace jtk;
  test_icp_1();
  test_icp_2();
  }