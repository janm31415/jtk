/*
   Do this:
      #define JTK_ICP_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_ICP_IMPLEMENTATION
   #include "jtk/icp.h"
 */

#ifndef JTK_ICP_H
#define JTK_ICP_H

#include "fitting.h"
#include "point_tree.h"
#include "mat.h"
#include "vec.h"

#ifndef JTKICPDEF
#ifdef JTK_ICP_STATIC
#define JTKICPDEF static
#else
#define JTKICPDEF extern
#endif
#endif

namespace jtk
  {

  JTKICPDEF matf16 icp(float& residual, 
    const std::vector<vec3<float>>& model_points,
    const std::vector<vec3<float>>& template_points,
    const matf16& initial_transformation,
    float inlier_distance=-1.f,
    uint32_t max_iterations=200,
    float tolerance=1e-4f,
    float inlier_shrink_factor=0.9f,
    float minimum_inlier_distance=0.05f,
    bool reflected = false
    );

  JTKICPDEF mat16 icp(double& residual,
    const std::vector<vec3<double>>& model_points,
    const std::vector<vec3<double>>& template_points,
    const mat16& initial_transformation,
    double inlier_distance = -1,
    uint32_t max_iterations = 200,
    double tolerance = 1e-4,
    double inlier_shrink_factor = 0.9,
    double minimum_inlier_distance = 0.05,
    bool reflected = false
  );

  } // namespace jtk

#endif //#ifndef JTK_ICP_H

#ifdef JTK_ICP_IMPLEMENTATION

#include <numeric>

namespace jtk
  {

  namespace
    {

    template <class T>
    struct point
      {
      vec3<T> pt;
      T& operator [] (size_t i)
        {
        return pt[i];
        }
      T operator [] (size_t i) const
        {
        return pt[i];
        }
      uint32_t idx;
      };

    template <class T>
    struct point_tree_traits
      {
      typedef T value_type;
      enum { dimension = 3 };
      typedef point<T> point_type;
      };

    template <class T>
    matrix<T, std::array<T, 16>> _icp(T& residual,
      const std::vector<vec3<T>>& model_points,
      const std::vector<vec3<T>>& template_points,
      const matrix<T, std::array<T, 16>>& initial_transformation,
      T inlier_distance,
      uint32_t max_iterations,
      T tolerance,
      T inlier_shrink_factor,
      T minimum_inlier_distance,
      bool reflected)
      {
      point_tree<point_tree_traits<T>> tree;
      std::vector<point<T>> pts;
      pts.reserve(model_points.size());
      for (uint32_t i = 0; i < (uint32_t)model_points.size(); ++i)
        {
        point<T> p;
        p.pt = model_points[i];
        p.idx = i;
        pts.push_back(p);
        }
      tree.efficient_build_tree(pts.begin(), pts.end());
      std::vector<uint32_t> active_template_points;
      if (inlier_distance <= 0)
        {
        active_template_points.resize(template_points.size());
        std::iota(active_template_points.begin(), active_template_points.end(), 0);
        }
      T delta = 1000;
      matrix<T, std::array<T, 16>> current = initial_transformation;
      for (uint32_t iter = 0; iter < max_iterations && delta > tolerance; ++iter)
        {
        if (inlier_distance > 0)
          {
          inlier_distance = std::max<T>(inlier_distance * inlier_shrink_factor, minimum_inlier_distance);
          active_template_points.clear();
          for (uint32_t i = 0; i < (uint32_t)template_points.size(); ++i)
            {
            matrix<T, std::array<T, 4>> v(4);
            v << template_points[i][0], template_points[i][1], template_points[i][2], 1;
            matrix<T, std::array<T, 4>> v_transf = current * v;
            point<T> vt;
            vt.pt = vec3<T>(v_transf(0), v_transf(1), v_transf(2));
            T dist;
            tree.find_nearest(dist, vt);
            if (dist < inlier_distance)
              active_template_points.push_back(i);
            }
          }

        matrix<T> source((uint32_t)active_template_points.size(), 3);
        matrix<T> destination((uint32_t)active_template_points.size(), 3);
        for (uint32_t i = 0; i < (uint32_t)active_template_points.size(); ++i)
          {
          matrix<T, std::array<T, 4>> v(4);
          v << template_points[active_template_points[i]][0], template_points[active_template_points[i]][1], template_points[active_template_points[i]][2], 1;
          matrix<T, std::array<T, 4>> v_transf = current * v;
          point<T> vt;
          vt.pt = vec3<T>(v_transf(0), v_transf(1), v_transf(2));
          T dist;
          auto nearest_source = tree.find_nearest(dist, vt);
          source(i, 0) = nearest_source.pt[0];
          source(i, 1) = nearest_source.pt[1];
          source(i, 2) = nearest_source.pt[2];
          destination(i, 0) = v_transf(0);
          destination(i, 1) = v_transf(1);
          destination(i, 2) = v_transf(2);
          }
        matrix<T, std::array<T, 16>> delta_transformation = reflected ? npoint_reflected(destination, source, false, true) : npoint(destination, source, false, true);
        current = delta_transformation * current;
        T l2_norm_rotation = (T)norm(block(delta_transformation, 0, 0, 3, 3) - identity(3, 3));
        T l2_norm_translation = (T)norm(block(delta_transformation, 0, 3, 3, 1));
        delta = std::max<T>(l2_norm_rotation, l2_norm_translation);
        }
      residual = 0.f;
      for (uint32_t i = 0; i < (uint32_t)active_template_points.size(); ++i)
        {
        matrix<T, std::array<T, 4>> v(4);
        v << template_points[active_template_points[i]][0], template_points[active_template_points[i]][1], template_points[active_template_points[i]][2], 1;
        matrix<T, std::array<T, 4>> v_transf = current * v;
        point<T> vt;
        vt.pt = vec3<T>(v_transf(0), v_transf(1), v_transf(2));
        T dist;
        auto nearest_source = tree.find_nearest(dist, vt);
        residual += dist;
        }
      residual /= (T)active_template_points.size();
      return current;
      }
    }

  JTKICPDEF matf16 icp(float& residual,
    const std::vector<vec3<float>>& model_points,
    const std::vector<vec3<float>>& template_points,
    const matf16& initial_transformation,
    float inlier_distance,
    uint32_t max_iterations,
    float tolerance,
    float inlier_shrink_factor,
    float minimum_inlier_distance,
    bool reflected)
    {
    return _icp<float>(residual, model_points, template_points, initial_transformation, inlier_distance, max_iterations, tolerance, inlier_shrink_factor, minimum_inlier_distance, reflected);    
    }

  JTKICPDEF mat16 icp(double& residual,
    const std::vector<vec3<double>>& model_points,
    const std::vector<vec3<double>>& template_points,
    const mat16& initial_transformation,
    double inlier_distance,
    uint32_t max_iterations,
    double tolerance,
    double inlier_shrink_factor,
    double minimum_inlier_distance,
    bool reflected)
    {
    return _icp<double>(residual, model_points, template_points, initial_transformation, inlier_distance, max_iterations, tolerance, inlier_shrink_factor, minimum_inlier_distance, reflected);
    }
  } // namespace jtk

#endif //JTK_ICP_IMPLEMENTATION