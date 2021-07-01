#pragma once

#include "mat.h"
#include "vec.h"

#include <stdexcept>

namespace jtk
  {

  ///////////////////////////////////////////////////////////////////////////////
  // interface
  ///////////////////////////////////////////////////////////////////////////////


  template <class T>
  matrix<T, std::array<T, 16>> npoint(const matrix<T>& source, const matrix<T>& destination, bool isotropic_scaling = false, bool correct_reflections = false);

  template <class T>
  matrix<T, std::array<T, 16>> npoint_negative(const matrix<T>& source, const matrix<T>& destination, bool isotropic_scaling = false, bool correct_reflections = false);

  template <class iterator>
  typename iterator::value_type centroid(iterator first, iterator last);

  template <class T>
  void fit_plane(jtk::vec3<T>& origin, jtk::vec3<T>& normal, T& smallest_eigenvalue, const std::vector<jtk::vec3<T>>& pts);

  ///////////////////////////////////////////////////////////////////////////////
  // implementation
  ///////////////////////////////////////////////////////////////////////////////

  namespace fitting_details
    {
    template <class T>
    matrix<T, std::array<T, 3>> mean(const matrix<T>& m)
      {
      matrix<T, std::array<T, 3>> out = zeros<T>(1, 3);
      for (uint32_t i = 0; i < m.rows(); ++i)
        {
        out(0, 0) += m(i, 0);
        out(0, 1) += m(i, 1);
        out(0, 2) += m(i, 2);
        }
      out(0, 0) /= (T)m.rows();
      out(0, 1) /= (T)m.rows();
      out(0, 2) /= (T)m.rows();
      return out;
      }
    } // namespace fitting_details

  template <class T>
  matrix<T, std::array<T, 16>> npoint(const matrix<T>& source, const matrix<T>& destination, bool isotropic_scaling, bool correct_reflections)
    {
    using namespace fitting_details;

    if (source.cols() != 3)
      throw std::runtime_error("npoint: source matrix should have 3 columns");
    if (destination.cols() != 3)
      throw std::runtime_error("npoint: destination matrix should have 3 columns");
    if (destination.rows() != source.rows())
      throw std::runtime_error("npoint: source and destination matrix should have equal amount of rows");

    matrix<T> X = destination;
    matrix<T> Y = source;
    auto Xmean = mean(destination);
    auto Ymean = mean(source);

    for (uint32_t i = 0; i < source.rows(); ++i)
      {
      X(i, 0) -= Xmean(0, 0);
      X(i, 1) -= Xmean(0, 1);
      X(i, 2) -= Xmean(0, 2);
      Y(i, 0) -= Ymean(0, 0);
      Y(i, 1) -= Ymean(0, 1);
      Y(i, 2) -= Ymean(0, 2);
      }

    matrix<T, std::array<T, 9>> A, U, VT;
    matrix<T, std::array<T, 3>> S;
    U = transpose(Y)*X;
    svd(U, S, VT);
    VT = transpose(VT);
    A = U * VT;

    if (correct_reflections && (determinant(A) < 0)) // special reflection case
      {
      VT(2, 0) *= -1.0;
      VT(2, 1) *= -1.0;
      VT(2, 2) *= -1.0;
      A = U * VT;
      }
    if (isotropic_scaling)
      A = U * VT*(S(0) + S(1) + S(2)) / pow(norm(Y), (T)2);

    matrix<T, std::array<T, 3>> trans = Xmean - Ymean * A;

    matrix<T, std::array<T, 16>> out(4, 4);

    for (uint32_t i = 0; i < 3; ++i)
      {
      for (uint32_t j = 0; j < 3; ++j)
        out(i, j) = A(j, i);
      out(i, 3) = trans(0, i);
      }
    out(3, 0) = (T)0;
    out(3, 1) = (T)0;
    out(3, 2) = (T)0;
    out(3, 3) = (T)1;

    return out;
    }

  template <class T>
  matrix<T, std::array<T, 16>> npoint_negative(const matrix<T>& source, const matrix<T>& destination, bool isotropic_scaling, bool correct_reflections)
    {
    using namespace fitting_details;

    if (source.cols() != 3)
      throw std::runtime_error("npoint: source matrix should have 3 columns");
    if (destination.cols() != 3)
      throw std::runtime_error("npoint: destination matrix should have 3 columns");
    if (destination.rows() != source.rows())
      throw std::runtime_error("npoint: source and destination matrix should have equal amount of rows");

    matrix<T> X = destination;
    matrix<T> Y = source;
    auto Xmean = mean(destination);
    auto Ymean = mean(source);

    for (uint32_t i = 0; i < source.rows(); ++i)
      {
      X(i, 0) -= Xmean(0, 0);
      X(i, 1) -= Xmean(0, 1);
      X(i, 2) -= Xmean(0, 2);
      Y(i, 0) -= Ymean(0, 0);
      Y(i, 1) -= Ymean(0, 1);
      Y(i, 2) -= Ymean(0, 2);
      }

    matrix<T, std::array<T, 9>> A, U, VT;
    matrix<T, std::array<T, 3>> S;
    U = transpose(Y) * X;
    svd(U, S, VT);
    VT = transpose(VT);
    VT(2, 0) *= -1.0;
    VT(2, 1) *= -1.0;
    VT(2, 2) *= -1.0;
    A = U * VT;

    if (correct_reflections && (determinant(A) > 0)) // special reflection case
      {
      VT(2, 0) *= -1.0;
      VT(2, 1) *= -1.0;
      VT(2, 2) *= -1.0;
      A = U * VT;
      }
    if (isotropic_scaling)
      A = U * VT * (S(0) + S(1) + S(2)) / pow(norm(Y), (T)2);

    matrix<T, std::array<T, 3>> trans = Xmean - Ymean * A;

    matrix<T, std::array<T, 16>> out(4, 4);

    for (size_t i = 0; i < 3; ++i)
      {
      for (size_t j = 0; j < 3; ++j)
        out(i, j) = A(j, i);
      out(i, 3) = trans(0, i);
      }
    out(3, 0) = (T)0;
    out(3, 1) = (T)0;
    out(3, 2) = (T)0;
    out(3, 3) = (T)1;

    return out;
    }

  template <class iterator>
  typename iterator::value_type centroid(iterator first, iterator last)
    {
    auto denom = std::distance(first, last);
    typename iterator::value_type pt((*first)[0], (*first)[1], (*first)[2]);
    ++first;
    for (; first != last; ++first)
      {
      pt[0] += ((*first)[0]);
      pt[1] += ((*first)[1]);
      pt[2] += ((*first)[2]);
      }
    pt[0] /= denom;
    pt[1] /= denom;
    pt[2] /= denom;
    return pt;
    }

  template <class T>
  void fit_plane(jtk::vec3<T>& origin, jtk::vec3<T>& normal, T& smallest_eigenvalue, const std::vector<jtk::vec3<T>>& pts)
    {
    origin = centroid(pts.begin(), pts.end());
    jtk::matrix<T, std::array<T, 9>> eigenvectors = jtk::zeros(3, 3);
    for (size_t k = 0; k < pts.size(); ++k)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          eigenvectors[i][j] += (pts[k][i] - origin[i])*(pts[k][j] - origin[j]);
    jtk::matrix<T, std::array<T, 9>> v(3, 3);
    jtk::matrix<T, std::array<T, 9>> eigenvalues(3, 1);
    svd(eigenvectors, eigenvalues, v);
    int smallest_eig = 0;
    if (std::abs(eigenvalues(1)) < std::abs(eigenvalues(smallest_eig)))
      smallest_eig = 1;
    if (std::abs(eigenvalues(2)) < std::abs(eigenvalues(smallest_eig)))
      smallest_eig = 2;
    smallest_eigenvalue = eigenvalues(smallest_eig);
    normal = jtk::vec3<T>(eigenvectors[0][smallest_eig], eigenvectors[1][smallest_eig], eigenvectors[2][smallest_eig]);
    }

  } // namespace jtk