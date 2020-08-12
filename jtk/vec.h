#pragma once

#include <cmath>

#undef min
#undef max

namespace jtk
  {
  template <typename T>
  struct vec2
    {
    T x, y;

    typedef T value_type;

    vec2();

    vec2(const T& t);

    vec2(const T& _x, const T& _y);

    template <typename T2>
    T& operator [] (T2 i);

    template <typename T2>
    const T& operator [] (T2 i) const;

    };

  template <typename T>
  struct vec3
    {
    T x, y, z;

    typedef T value_type;

    vec3();
    vec3(const T& t);

    vec3(const T& _x, const T& _y, const T& _z);

    template <typename T2>
    T& operator [] (T2 i);

    template <typename T2>
    const T& operator [] (T2 i) const;

    };


  template <typename T>
  struct vec4
    {
    T x, y, z, w;

    typedef T value_type;

    vec4();

    vec4(const T& t);

    vec4(const T& _x, const T& _y, const T& _z, const T& _w);

    template <typename T2>
    T& operator [] (T2 i);

    template <typename T2>
    const T& operator [] (T2 i) const;

    };

  template <class T>
  struct boundingbox3d
    {
    vec3<T> min, max;
    };

  template <class T, class Iter>
  boundingbox3d<T> bounding_volume_3d(Iter first, Iter last);

  template <class T>
  inline bool inside(const boundingbox3d<T>& bb, const vec3<T>& pt);

  template <class T>
  inline void add_point(boundingbox3d<T>& bb, const vec3<T>& pt);

  template <class T>
  inline vec3<T> center(const boundingbox3d<T>& bb);

  template <class T>
  inline boundingbox3d<T> unite(const boundingbox3d<T>& left, const boundingbox3d<T>& right);

  template <class T>
  struct boundingbox2d
    {
    vec2<T> min, max;
    };

  template <class T, class Iter>
  boundingbox2d<T> bounding_volume_2d(Iter first, Iter last);

  template <class T>
  inline bool inside(const boundingbox2d<T>& bb, const vec2<T>& pt);

  template <class T>
  inline void add_point(boundingbox2d<T>& bb, const vec2<T>& pt);

  template <class T>
  inline vec2<T> center(const boundingbox2d<T>& bb);

  template <class T>
  inline boundingbox2d<T> unite(const boundingbox2d<T>& left, const boundingbox2d<T>& right);

  /////////////////////////////////////////////////////////////////////////
  // struct vec2
  /////////////////////////////////////////////////////////////////////////

  template <typename T>
  vec2<T>::vec2() {}

  template <typename T>
  vec2<T>::vec2(const T& t) : x(t), y(t) {}

  template <typename T>
  vec2<T>::vec2(const T& _x, const T& _y) : x(_x), y(_y) {}

  template <typename T>
  template <typename T2>
  T& vec2<T>::operator [] (T2 i) { return (&x)[i]; }

  template <typename T>
  template <typename T2>
  const T& vec2<T>::operator [] (T2 i) const { return (&x)[i]; }

  template <typename T>
  inline vec2<T> operator + (const vec2<T>& a)
    {
    return vec2<T>(+a.x, +a.y);
    }

  template <typename T>
  inline vec2<T> operator - (const vec2<T>& a)
    {
    return vec2<T>(-a.x, -a.y);
    }

  template <typename T>
  inline vec2<T> operator + (const vec2<T>& left, const vec2<T>& right)
    {
    return vec2<T>(left.x + right.x, left.y + right.y);
    }

  template <typename T>
  inline vec2<T> operator - (const vec2<T>& left, const vec2<T>& right)
    {
    return vec2<T>(left.x - right.x, left.y - right.y);
    }

  template <typename T>
  inline vec2<T> operator * (const vec2<T>& left, const vec2<T>& right)
    {
    return vec2<T>(left.x * right.x, left.y * right.y);
    }

  template <typename T>
  inline vec2<T> operator * (const T& s, const vec2<T>& right)
    {
    return vec2<T>(s * right.x, s * right.y);
    }

  template <typename T>
  inline vec2<T> operator * (const vec2<T>& left, const T& s)
    {
    return vec2<T>(left.x * s, left.y * s);
    }

  template <typename T>
  inline vec2<T> operator / (const vec2<T>& left, const vec2<T>& right)
    {
    return vec2<T>(left.x / right.x, left.y / right.y);
    }

  template <typename T>
  inline vec2<T> operator / (const T& s, const vec2<T>& right)
    {
    return vec2<T>(s / right.x, s / right.y);
    }

  template <typename T>
  inline vec2<T> operator / (const vec2<T>& left, const T& s)
    {
    return vec2<T>(left.x / s, left.y / s);
    }

  template <typename T>
  inline vec2<T> min(const vec2<T>& left, const vec2<T>& right)
    {
    return vec2<T>(min(left.x, right.x), min(left.y, right.y));
    }

  template <typename T>
  inline vec2<T> max(const vec2<T>& left, const vec2<T>& right)
    {
    return vec2<T>(max(left.x, right.x), max(left.y, right.y));
    }

  template <typename T>
  inline vec2<T> abs(const vec2<T>& a)
    {
    return vec2<T>(abs(a.x), abs(a.y));
    }

  template <>
  inline vec2<float> abs(const vec2<float>& a)
    {
    return vec2<float>(fabs(a.x), fabs(a.y));
    }

  template <typename T>
  inline vec2<T> sqrt(const vec2<T>& a)
    {
    return vec2<T>(sqrt(a.x), sqrt(a.y));
    }

  template <typename T>
  inline vec2<T> rsqrt(const vec2<T>& a)
    {
    return vec2<T>(rsqrt(a.x), rsqrt(a.y));
    }

  template <typename T>
  inline vec2<T> reciprocal(const vec2<T>& a)
    {
    return vec2<T>(reciprocal(a.x), reciprocal(a.y));
    }

  template <typename T>
  inline T dot(const vec2<T>& left, const vec2<T>& right)
    {
    return left.x*right.x + left.y*right.y;
    }

  template <typename T>
  inline T length(const vec2<T>& a)
    {
    return sqrt(dot(a, a));
    }

  template <>
  inline float length(const vec2<float>& a)
    {
    return std::sqrt(dot(a, a));
    }

  template <>
  inline double length(const vec2<double>& a)
    {
    return std::sqrt(dot(a, a));
    }

  template <typename T>
  inline T distance(const vec2<T>& left, const vec2<T>& right)
    {
    return length(left - right);
    }

  template <typename T>
  inline T distance_sqr(const vec2<T>& left, const vec2<T>& right)
    {
    auto a = (left - right);
    return dot(a, a);
    }

  template <typename T>
  inline vec2<T> normalize(const vec2<T>& a)
    {
    return a * rsqrt(dot(a, a));
    }

  template <>
  inline vec2<float> normalize(const vec2<float>& a)
    {
    return a / std::sqrt(dot(a, a));
    }

  template <>
  inline vec2<double> normalize(const vec2<double>& a)
    {
    return a / std::sqrt(dot(a, a));
    }

  template <typename T>
  inline bool operator < (const vec2<T>& left, const vec2<T>& right)
    {
    return (left[0] == right[0]) ? (left[1] < right[1]) : (left[0] < right[0]);
    }

  template <typename T>
  inline bool operator > (const vec2<T>& left, const vec2<T>& right)
    {
    return (left[0] == right[0]) ? (left[1] > right[1]) : (left[0] > right[0]);
    }

  template <typename T>
  inline bool operator == (const vec2<T>& left, const vec2<T>& right)
    {
    return (left[0] == right[0]) && (left[1] == right[1]);
    }

  template <typename T>
  inline bool operator != (const vec2<T>& left, const vec2<T>& right)
    {
    return !(left == right);
    }

  /////////////////////////////////////////////////////////////////////////
  // struct vec3
  /////////////////////////////////////////////////////////////////////////

  template <typename T>
  vec3<T>::vec3() {}

  template <typename T>
  vec3<T>::vec3(const T& t) : x(t), y(t), z(t) {}

  template <typename T>
  vec3<T>::vec3(const T& _x, const T& _y, const T& _z) : x(_x), y(_y), z(_z) {}

  template <typename T>
  template <typename T2>
  T& vec3<T>::operator [] (T2 i) { return (&x)[i]; }

  template <typename T>
  template <typename T2>
  const T& vec3<T>::operator [] (T2 i) const { return (&x)[i]; }

  template <typename T>
  inline vec3<T> operator + (const vec3<T>& a)
    {
    return vec3<T>(+a.x, +a.y, +a.z);
    }

  template <typename T>
  inline vec3<T> operator - (const vec3<T>& a)
    {
    return vec3<T>(-a.x, -a.y, -a.z);
    }

  template <typename T>
  inline vec3<T> operator + (const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(left.x + right.x, left.y + right.y, left.z + right.z);
    }

  template <typename T>
  inline vec3<T> operator - (const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(left.x - right.x, left.y - right.y, left.z - right.z);
    }

  template <typename T>
  inline vec3<T> operator * (const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(left.x * right.x, left.y * right.y, left.z * right.z);
    }

  template <typename T>
  inline vec3<T> operator * (const T& s, const vec3<T>& right)
    {
    return vec3<T>(s * right.x, s * right.y, s * right.z);
    }

  template <typename T>
  inline vec3<T> operator * (const vec3<T>& left, const T& s)
    {
    return vec3<T>(left.x * s, left.y * s, left.z * s);
    }

  template <typename T>
  inline vec3<T> operator / (const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(left.x / right.x, left.y / right.y, left.z / right.z);
    }

  template <typename T>
  inline vec3<T> operator / (const T& s, const vec3<T>& right)
    {
    return vec3<T>(s / right.x, s / right.y, s / right.z);
    }

  template <typename T>
  inline vec3<T> operator / (const vec3<T>& left, const T& s)
    {
    return vec3<T>(left.x / s, left.y / s, left.z / s);
    }

  template <typename T>
  inline vec3<T> min(const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(min(left.x, right.x), min(left.y, right.y), min(left.z, right.z));
    }

  template <typename T>
  inline vec3<T> max(const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(max(left.x, right.x), max(left.y, right.y), max(left.z, right.z));
    }

  template <>
  inline vec3<double> min(const vec3<double>& left, const vec3<double>& right)
    {
    return vec3<double>(left.x < right.x ? left.x : right.x, left.y < right.y ? left.y : right.y, left.z < right.z ? left.z : right.z);
    }

  template <>
  inline vec3<double> max(const vec3<double>& left, const vec3<double>& right)
    {
    return vec3<double>(left.x < right.x ? right.x : left.x, left.y < right.y ? right.y : left.y, left.z < right.z ? right.z : left.z);
    }

  template <>
  inline vec3<float> min(const vec3<float>& left, const vec3<float>& right)
    {
    return vec3<float>(left.x < right.x ? left.x : right.x, left.y < right.y ? left.y : right.y, left.z < right.z ? left.z : right.z);
    }

  template <>
  inline vec3<float> max(const vec3<float>& left, const vec3<float>& right)
    {
    return vec3<float>(left.x < right.x ? right.x : left.x, left.y < right.y ? right.y : left.y, left.z < right.z ? right.z : left.z);
    }

  template <typename T>
  inline vec3<T> abs(const vec3<T>& a)
    {
    return vec3<T>(abs(a.x), abs(a.y), abs(a.z));
    }

  template <>
  inline vec3<float> abs(const vec3<float>& a)
    {
    return vec3<float>(fabs(a.x), fabs(a.y), fabs(a.z));
    }

  template <typename T>
  inline vec3<T> sqrt(const vec3<T>& a)
    {
    return vec3<T>(sqrt(a.x), sqrt(a.y), sqrt(a.z));
    }

  template <>
  inline vec3<float> sqrt(const vec3<float>& a)
    {
    return vec3<float>(std::sqrt(a.x), std::sqrt(a.y), std::sqrt(a.z));
    }

  template <typename T>
  inline vec3<T> rsqrt(const vec3<T>& a)
    {
    return vec3<T>(rsqrt(a.x), rsqrt(a.y), rsqrt(a.z));
    }

  template <typename T>
  inline vec3<T> reciprocal(const vec3<T>& a)
    {
    return vec3<T>(reciprocal(a.x), reciprocal(a.y), reciprocal(a.z));
    }

  template <>
  inline vec3<float> reciprocal<float>(const vec3<float>& a)
    {
    return vec3<float>(1.f / a.x, 1.f / a.y, 1.f / a.z);
    }

  template <typename T>
  inline vec3<T> cross(const vec3<T>& left, const vec3<T>& right)
    {
    return vec3<T>(left.y*right.z - left.z*right.y, left.z*right.x - left.x*right.z, left.x*right.y - left.y*right.x);
    }

  template <typename T>
  inline T dot(const vec3<T>& left, const vec3<T>& right)
    {
    return left.x*right.x + left.y*right.y + left.z*right.z;
    }

  template <typename T>
  inline T length_sqr(const vec3<T>& a)
    {
    return dot(a, a);
    }

  template <typename T>
  inline T length(const vec3<T>& a)
    {
    return sqrt(dot(a, a));
    }

  template <>
  inline float length(const vec3<float>& a)
    {
    return std::sqrt(dot(a, a));
    }

  template <>
  inline double length(const vec3<double>& a)
    {
    return std::sqrt(dot(a, a));
    }

  template <typename T>
  inline T distance(const vec3<T>& left, const vec3<T>& right)
    {
    return length(left - right);
    }

  template <typename T>
  inline T distance_sqr(const vec3<T>& left, const vec3<T>& right)
    {
    return length_sqr(left - right);
    }

  template <typename T>
  inline vec3<T> normalize(const vec3<T>& a)
    {
    return a * rsqrt(dot(a, a));
    }

  template <>
  inline vec3<float> normalize(const vec3<float>& a)
    {
    float denom = std::sqrt(dot(a, a));
    return denom ? a / denom : a;
    }

  template <>
  inline vec3<double> normalize(const vec3<double>& a)
    {
    double denom = std::sqrt(dot(a, a));
    return denom ? a / denom : a;
    }

  template <typename T>
  inline bool operator < (const vec3<T>& left, const vec3<T>& right)
    {
    return (left[0] == right[0]) ? ((left[1] == right[1]) ? (left[2] < right[2]) : (left[1] < right[1])) : (left[0] < right[0]);
    }

  template <typename T>
  inline bool operator > (const vec3<T>& left, const vec3<T>& right)
    {
    return (left[0] == right[0]) ? ((left[1] == right[1]) ? (left[2] > right[2]) : (left[1] > right[1])) : (left[0] > right[0]);
    }

  template <typename T>
  inline bool operator == (const vec3<T>& left, const vec3<T>& right)
    {
    return (left[0] == right[0]) && (left[1] == right[1]) && (left[2] == right[2]);
    }

  template <typename T>
  inline bool operator != (const vec3<T>& left, const vec3<T>& right)
    {
    return !(left == right);
    }

  /////////////////////////////////////////////////////////////////////////
  // struct vec4
  /////////////////////////////////////////////////////////////////////////

  template <typename T>
  vec4<T>::vec4() {}

  template <typename T>
  vec4<T>::vec4(const T& t) : x(t), y(t), z(t), w(t) {}

  template <typename T>
  vec4<T>::vec4(const T& _x, const T& _y, const T& _z, const T& _w) : x(_x), y(_y), z(_z), w(_w) {}

  template <typename T>
  template <typename T2>
  T& vec4<T>::operator [] (T2 i) { return (&x)[i]; }

  template <typename T>
  template <typename T2>
  const T& vec4<T>::operator [] (T2 i) const { return (&x)[i]; }

  template <typename T>
  inline vec4<T> operator + (const vec4<T>& a)
    {
    return vec4<T>(+a.x, +a.y, +a.z, +a.w);
    }

  template <typename T>
  inline vec4<T> operator - (const vec4<T>& a)
    {
    return vec4<T>(-a.x, -a.y, -a.z, -a.w);
    }

  template <typename T>
  inline vec4<T> operator + (const vec4<T>& left, const vec4<T>& right)
    {
    return vec4<T>(left.x + right.x, left.y + right.y, left.z + right.z, left.w + right.w);
    }

  template <typename T>
  inline vec4<T> operator - (const vec4<T>& left, const vec4<T>& right)
    {
    return vec4<T>(left.x - right.x, left.y - right.y, left.z - right.z, left.w - right.w);
    }

  template <typename T>
  inline vec4<T> operator * (const vec4<T>& left, const vec4<T>& right)
    {
    return vec4<T>(left.x * right.x, left.y * right.y, left.z * right.z, left.w * right.w);
    }

  template <typename T>
  inline vec4<T> operator * (const T& s, const vec4<T>& right)
    {
    return vec4<T>(s * right.x, s * right.y, s * right.z, s * right.w);
    }

  template <typename T>
  inline vec4<T> operator * (const vec4<T>& left, const T& s)
    {
    return vec4<T>(left.x * s, left.y * s, left.z * s, left.w * s);
    }

  template <typename T>
  inline vec4<T> operator / (const vec4<T>& left, const vec4<T>& right)
    {
    return vec4<T>(left.x / right.x, left.y / right.y, left.z / right.z, left.w / right.w);
    }

  template <typename T>
  inline vec4<T> operator / (const T& s, const vec4<T>& right)
    {
    return vec4<T>(s / right.x, s / right.y, s / right.z, s / right.w);
    }

  template <typename T>
  inline vec4<T> operator / (const vec4<T>& left, const T& s)
    {
    return vec4<T>(left.x / s, left.y / s, left.z / s, left.w / s);
    }

  template <typename T>
  inline vec4<T> min(const vec4<T>& left, const vec4<T>& right)
    {
    return vec4<T>(min(left.x, right.x), min(left.y, right.y), min(left.z, right.z), min(left.w, right.w));
    }

  template <typename T>
  inline vec4<T> max(const vec4<T>& left, const vec4<T>& right)
    {
    return vec4<T>(max(left.x, right.x), max(left.y, right.y), max(left.z, right.z), max(left.w, right.w));
    }

  template <>
  inline vec4<double> min(const vec4<double>& left, const vec4<double>& right)
    {
    return vec4<double>(left.x < right.x ? left.x : right.x, left.y < right.y ? left.y : right.y, left.z < right.z ? left.z : right.z, left.w < right.w ? left.w : right.w);
    }

  template <>
  inline vec4<double> max(const vec4<double>& left, const vec4<double>& right)
    {
    return vec4<double>(left.x < right.x ? right.x : left.x, left.y < right.y ? right.y : left.y, left.z < right.z ? right.z : left.z, left.w < right.w ? right.w : left.w);
    }

  template <>
  inline vec4<float> min(const vec4<float>& left, const vec4<float>& right)
    {
    return vec4<float>(left.x < right.x ? left.x : right.x, left.y < right.y ? left.y : right.y, left.z < right.z ? left.z : right.z, left.w < right.w ? left.w : right.w);
    }

  template <>
  inline vec4<float> max(const vec4<float>& left, const vec4<float>& right)
    {
    return vec4<float>(left.x < right.x ? right.x : left.x, left.y < right.y ? right.y : left.y, left.z < right.z ? right.z : left.z, left.w < right.w ? right.w : left.w);
    }

  template <typename T>
  inline vec4<T> abs(const vec4<T>& a)
    {
    return vec4<T>(abs(a.x), abs(a.y), abs(a.z), abs(a.w));
    }

  template <>
  inline vec4<float> abs(const vec4<float>& a)
    {
    return vec4<float>(fabs(a.x), fabs(a.y), fabs(a.z), fabs(a.w));
    }

  template <typename T>
  inline vec4<T> sqrt(const vec4<T>& a)
    {
    return vec4<T>(sqrt(a.x), sqrt(a.y), sqrt(a.z), sqrt(a.w));
    }

  template <typename T>
  inline vec4<T> rsqrt(const vec4<T>& a)
    {
    return vec4<T>(rsqrt(a.x), rsqrt(a.y), rsqrt(a.z), rsqrt(a.w));
    }

  template <typename T>
  inline vec4<T> reciprocal(const vec4<T>& a)
    {
    return vec4<T>(reciprocal(a.x), reciprocal(a.y), reciprocal(a.z), reciprocal(a.w));
    }

  template <>
  inline vec4<float> reciprocal<float>(const vec4<float>& a)
    {
    return vec4<float>(1.f / a.x, 1.f / a.y, 1.f / a.z, 1.f / a.w);
    }

  template <typename T>
  inline bool operator == (const vec4<T>& left, const vec4<T>& right)
    {
    return (left[0] == right[0]) && (left[1] == right[1]) && (left[2] == right[2]) && (left[3] == right[3]);
    }

  template <typename T>
  inline bool operator != (const vec4<T>& left, const vec4<T>& right)
    {
    return !(left == right);
    }

  /////////////////////////////////////////////////////////////////////////
  // bounding 
  /////////////////////////////////////////////////////////////////////////

  template <class T, class Iter>
  boundingbox3d<T> bounding_volume_3d(Iter first, Iter last)
    {
    boundingbox3d<T> result;
    if (first == last)
      {
      result.min[0] = (T)1;
      result.min[1] = (T)1;
      result.min[2] = (T)1;
      result.max[0] = (T)0;
      result.max[1] = (T)0;
      result.max[2] = (T)0;
      return result;
      }
    result.min[0] = (T)((*first)[0]);
    result.min[1] = (T)((*first)[1]);
    result.min[2] = (T)((*first)[2]);
    result.max[0] = (T)((*first)[0]);
    result.max[1] = (T)((*first)[1]);
    result.max[2] = (T)((*first)[2]);
    ++first;
    for (; first != last; ++first)
      {
      const auto& pt = *first;
      for (int i = 0; i < 3; ++i)
        {
        if (result.min[i] > (T)(pt[i]))
          result.min[i] = (T)(pt[i]);
        if (result.max[i] < (T)(pt[i]))
          result.max[i] = (T)(pt[i]);
        }
      }
    return result;
    }

  template <class T>
  inline bool inside(const boundingbox3d<T>& bb, const vec3<T>& pt)
    {
    for (size_t i = 0; i < 3; ++i)
      {
      if (pt[i] < bb.min[i])
        return false;
      if (pt[i] > bb.max[i])
        return false;
      }
    return true;
    }

  template <class T>
  inline void add_point(boundingbox3d<T>& bb, const vec3<T>& pt)
    {
    for (size_t i = 0; i < 3; ++i)
      {
      if (bb.min[i] > pt[i])
        bb.min[i] = pt[i];
      if (bb.max[i] < pt[i])
        bb.max[i] = pt[i];
      }
    }

  template <class T>
  inline vec3<T> center(const boundingbox3d<T>& bb)
    {
    vec3<T> c;
    c[0] = (bb.min[0] + bb.max[0]) / 2;
    c[1] = (bb.min[1] + bb.max[1]) / 2;
    c[2] = (bb.min[2] + bb.max[2]) / 2;
    return c;
    }

  template <class T>
  inline boundingbox3d<T> unite(const boundingbox3d<T>& left, const boundingbox3d<T>& right)
    {
    boundingbox3d<T> bb(left);
    if (right.min[0] <= right.max[0] &&
      right.min[1] <= right.max[1] &&
      right.min[2] <= right.max[2])
      {
      add_point(bb, right.min);
      add_point(bb, right.max);
      }
    return bb;
    }

  template <class T, class Iter>
  boundingbox2d<T> bounding_volume_2d(Iter first, Iter last)
    {
    boundingbox2d<T> result;
    if (first == last)
      {
      result.min[0] = (T)1;
      result.min[1] = (T)1;
      result.max[0] = (T)0;
      result.max[1] = (T)0;
      return result;
      }
    result.min[0] = (T)((*first)[0]);
    result.min[1] = (T)((*first)[1]);
    result.max[0] = (T)((*first)[0]);
    result.max[1] = (T)((*first)[1]);
    ++first;
    for (; first != last; ++first)
      {
      const auto& pt = *first;
      for (int i = 0; i < 2; ++i)
        {
        if (result.min[i] > (T)(pt[i]))
          result.min[i] = (T)(pt[i]);
        if (result.max[i] < (T)(pt[i]))
          result.max[i] = (T)(pt[i]);
        }
      }
    return result;
    }

  template <class T>
  inline bool inside(const boundingbox2d<T>& bb, const vec2<T>& pt)
    {
    for (size_t i = 0; i < 2; ++i)
      {
      if (pt[i] < bb.min[i])
        return false;
      if (pt[i] > bb.max[i])
        return false;
      }
    return true;
    }

  template <class T>
  inline void add_point(boundingbox2d<T>& bb, const vec2<T>& pt)
    {
    for (size_t i = 0; i < 2; ++i)
      {
      if (bb.min[i] > pt[i])
        bb.min[i] = pt[i];
      if (bb.max[i] < pt[i])
        bb.max[i] = pt[i];
      }
    }

  template <class T>
  inline vec2<T> center(const boundingbox2d<T>& bb)
    {
    vec2<T> c;
    c[0] = (bb.min[0] + bb.max[0]) / 2;
    c[1] = (bb.min[1] + bb.max[1]) / 2;
    return c;
    }

  template <class T>
  inline boundingbox2d<T> unite(const boundingbox2d<T>& left, const boundingbox2d<T>& right)
    {
    boundingbox2d<T> bb(left);
    if (right.min[0] <= right.max[0] &&
      right.min[1] <= right.max[1])
      {
      add_point(bb, right.min);
      add_point(bb, right.max);
      }
    return bb;
    }

  } // namespace jtk