/*
   Do this:
      #define JTK_QBVH_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_QBVH_IMPLEMENTATION
   #include "jtk/qbvh.h"
 */

#ifndef JTK_QBVH_H
#define JTK_QBVH_H

#include "concurrency.h"
#include "vec.h"

#ifdef _JTK_FOR_ARM
#include "sse2neon.h"
#else
#include <immintrin.h>
#include <emmintrin.h>
#endif

#include <array>
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <string.h>

#ifndef JTKQBVHDEF
#ifdef JTK_QBVH_STATIC
#define JTKQBVHDEF static
#else
#define JTKQBVHDEF extern
#endif
#endif

#ifndef JTKQBVHINLINE
#define JTKQBVHINLINE inline
#endif

namespace jtk
  {
  /////////////////////////////////////////////////////////////////////////
  // interfaces
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF void* aligned_malloc(size_t size, size_t align);
  JTKQBVHDEF void aligned_free(void* ptr);
  template<typename T, size_t alignment>
  struct aligned_allocator
    {
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    pointer allocate(size_type n);
    void deallocate(pointer p, size_type);
    void construct(pointer p, const_reference val);
    void destroy(pointer p);
    };

  template<typename T, typename allocator>
  class vector_t
    {
    public:
      typedef T value_type;
      typedef T* iterator;
      typedef const T* const_iterator;

      vector_t();
      explicit vector_t(size_t sz);
      explicit vector_t(size_t sz, const allocator& alloc);
      vector_t(size_t sz, const T& val, const allocator& alloc = allocator());
      ~vector_t();
      vector_t(const vector_t& other);
      vector_t(vector_t&& other);
      vector_t& operator=(const vector_t& other);
      vector_t& operator=(vector_t&& other);
      void swap(vector_t& other);
            iterator begin();
      const_iterator begin() const;
            iterator end();
      const_iterator end() const;
      bool   empty() const;
      size_t size() const;
      size_t capacity() const;
      void resize(size_t new_size);
      void reserve(size_t new_alloced);
      void shrink_to_fit();
            T& operator[](size_t i);
      const T& operator[](size_t i) const;
            T& at(size_t i);
      const T& at(size_t i) const;
      T& front() const;
      T& back() const;
            T* data();
      const T* data() const;
      void push_back(const T& nt);
      void pop_back();
      void clear();
      friend bool operator== (const vector_t& a, const vector_t& b);
      friend bool operator!= (const vector_t& a, const vector_t& b);

    private:
      void internal_resize_init(size_t new_active);
      void internal_resize_fill(size_t new_active, const T& val);
      void internal_resize(size_t new_active, size_t new_alloced);
      size_t internal_grow_size(size_t new_alloced);

    private:
      allocator alloc;
      size_t size_active;    // number of valid items
      size_t size_alloced;   // number of items allocated
      T* items;              // data array
    };

  template<typename T>
  using aligned_vector = vector_t<T, aligned_allocator<T, std::alignment_of<T>::value> >;

#ifdef _WIN32
  __declspec(align(16))
#endif
    struct bool4
    {
    union
      {
      __m128 m128;
      uint32_t u[4];
      int32_t i[4];
      };

    bool4();
    bool4(const __m128 in);
    bool4(const __m128i in);
    bool4(bool b);
    bool4(bool b0, bool b1, bool b2, bool b3);

    template <typename T>
    int32_t& operator [](T index);
    template <typename T>
    int32_t operator [](T index) const;
    }
#ifndef _WIN32 // linux alignment in gcc
  __attribute__((aligned(16)))
#endif
    ;
  JTKQBVHDEF bool all(const bool4& b);
  JTKQBVHDEF bool any(const bool4& b);
  JTKQBVHDEF bool none(const bool4& b);
  JTKQBVHDEF bool4 operator ! (const bool4& a);
  JTKQBVHDEF bool4 operator & (const bool4& a, const bool4& b);
  JTKQBVHDEF bool4 operator | (const bool4& a, const bool4& b);
  JTKQBVHDEF bool4 operator ^ (const bool4& a, const bool4& b);
  JTKQBVHDEF bool4& operator &= (bool4& a, const bool4& b);
  JTKQBVHDEF bool4& operator |= (bool4& a, const bool4& b);
  JTKQBVHDEF bool4& operator ^= (bool4& a, const bool4& b);
  JTKQBVHDEF bool4 operator != (const bool4& a, const bool4& b);
  JTKQBVHDEF bool4 operator == (const bool4& a, const bool4& b);

#ifdef _WIN32
  __declspec(align(16))
#endif
    struct int4
    {
    union
      {
      __m128i m128i;
      int32_t i[4];
      uint32_t u[4];
      };

    template <typename T>
    int32_t& operator [] (T ind);

    template <typename T>
    int32_t operator [] (T ind) const;

    int4();
    int4(const __m128i& in);
    int4(int32_t in);
    int4(int32_t i0, int32_t i1, int32_t i2, int32_t i3);
    int4(const bool4& in);

    typedef int32_t value_type;
    }
#ifndef _WIN32 // linux alignment in gcc
  __attribute__((aligned(16)))
#endif  
    ;

  JTKQBVHDEF int4 operator + (const int4& a);
  JTKQBVHDEF int4 operator - (const int4& a);
  JTKQBVHDEF int4 operator + (const int4& left, const int4& right);
  JTKQBVHDEF int4 operator - (const int4& left, const int4& right);
  JTKQBVHDEF int4 operator * (const int4& left, const int4& right);
  JTKQBVHDEF int4 operator * (const int4& left, int32_t right);
  JTKQBVHDEF int4 operator * (int32_t left, const int4& right);
  JTKQBVHDEF int4 min(const int4& left, const int4& right);
  JTKQBVHDEF int4 max(const int4& left, const int4& right);
  JTKQBVHDEF bool4 operator == (const int4& left, const int4& right);
  JTKQBVHDEF bool4 operator != (const int4& left, const int4& right);
  JTKQBVHDEF bool4 operator < (const int4& left, const int4& right);
  JTKQBVHDEF bool4 operator > (const int4& left, const int4& right);
  JTKQBVHDEF bool4 operator <= (const int4& left, const int4& right);
  JTKQBVHDEF bool4 operator >= (const int4& left, const int4& right);
  JTKQBVHDEF int4 masked_update(const bool4& mask, const int4& original, const int4& updated_values);
  JTKQBVHDEF int4 operator & (const int4& left, const int4& right);
  JTKQBVHDEF int4 operator | (const int4& left, const int4& right);
#ifndef _JTK_FOR_ARM
  JTKQBVHDEF int4 operator >> (const int4& a, int n);
  JTKQBVHDEF int4 operator << (const int4& a, int n);
#endif
  JTKQBVHDEF int any(const int4& a);

#ifdef _WIN32
  __declspec(align(16))
#endif
    struct float4
    {
    union
      {
      __m128 m128;
      float f[4];
      uint32_t u[4];
      int32_t i[4];
      };

    template <typename T>
    float& operator [] (T i);
    template <typename T>
    float operator [] (T i) const;

    float4();
    float4(const __m128 in);
    float4(float f);
    float4(float _x, float _y, float _z);
    float4(float _x, float _y, float _z, float _w);

    typedef float value_type;
    }
#ifndef _WIN32 // linux alignment in gcc
  __attribute__((aligned(16)))
#endif  
    ;
  JTKQBVHDEF float4 operator + (const float4& a);
  JTKQBVHDEF float4 operator - (const float4& a);
  JTKQBVHDEF float4 operator + (const float4& left, const float4& right);
  JTKQBVHDEF float4 operator - (const float4& left, const float4& right);
  JTKQBVHDEF float4 operator * (const float4& left, const float4& right);
  JTKQBVHDEF float4 operator * (const float4& left, float right);
  JTKQBVHDEF float4 operator * (float left, const float4& right);
  JTKQBVHDEF float4 operator / (const float4& left, const float4& right);
  JTKQBVHDEF float4 operator / (const float4& left, float right);
  JTKQBVHDEF float4 operator / (float left, const float4& right);
  JTKQBVHDEF float4 min(const float4& left, const float4& right);
  JTKQBVHDEF float4 max(const float4& left, const float4& right);
  JTKQBVHDEF float min_horizontal(const float4& x);
  JTKQBVHDEF float max_horizontal(const float4& x);
  JTKQBVHDEF float4 cross(const float4& left, const float4& right);
  JTKQBVHDEF float dot(const float4& left, const float4& right);
  JTKQBVHDEF float dot4(const float4& left, const float4& right);
  JTKQBVHDEF float4 abs(const float4& a);
  JTKQBVHDEF float4 sqrt(const float4& a);
  JTKQBVHDEF float4 rsqrt(const float4& a);
  JTKQBVHDEF float4 reciprocal(const float4& a);
  JTKQBVHDEF bool4 operator == (const float4& left, const float4& right);
  JTKQBVHDEF bool4 operator != (const float4& left, const float4& right);
  JTKQBVHDEF bool4 operator < (const float4& left, const float4& right);
  JTKQBVHDEF bool4 operator > (const float4& left, const float4& right);
  JTKQBVHDEF bool4 operator <= (const float4& left, const float4& right);
  JTKQBVHDEF bool4 operator >= (const float4& left, const float4& right);
  JTKQBVHDEF float4 unpacklo(const float4& left, const float4& right);
  JTKQBVHDEF float4 unpackhi(const float4& left, const float4& right);
  JTKQBVHDEF void transpose(float4& r0, float4& r1, float4& r2, float4& r3, const float4& c0, const float4& c1, const float4& c2, const float4& c3);
  JTKQBVHDEF float4 masked_update(const bool4& mask, const float4& original, const float4& updated_values);
  JTKQBVHDEF float4 masked_update(const int4& mask, const float4& original, const float4& updated_values);


#ifdef _WIN32
  _declspec(align(16))
#endif
    struct float4x4
    {
    union
      {
      float4 col[4];
      float f[16];
      };

    template <typename T>
    float& operator [] (T i);

    template <typename T>
    float operator [] (T i) const;

    float4x4();
    float4x4(const float4& col0, const float4& col1, const float4& col2, const float4& col3);
    float4x4(float* m);
    }
#ifndef _WIN32 // linux alignment in gcc
  __attribute__((aligned(16)))
#endif  
    ;
  JTKQBVHDEF float4x4 get_identity();
  JTKQBVHDEF float4x4 make_translation(float x, float y, float z);
  JTKQBVHDEF float4x4 transpose(const float4x4& m);
  JTKQBVHDEF float4x4 invert_orthonormal(const float4x4& m);
  // for column major matrix
  // we use __m128 to represent 2x2 matrix as A = | A0  A1 |
  //                                              | A2  A3 |
  // 2x2 column major matrix multiply A*B
  JTKQBVHDEF __m128 mat2mul(__m128 vec1, __m128 vec2);
  // 2x2 column major matrix adjugate multiply (A#)*B
  JTKQBVHDEF __m128 mat2adjmul(__m128 vec1, __m128 vec2);
  // 2x2 column major matrix multiply adjugate A*(B#)
  JTKQBVHDEF __m128 mat2muladj(__m128 vec1, __m128 vec2);
  JTKQBVHDEF float4x4 invert(const float4x4& m);
  JTKQBVHDEF float4 matrix_vector_multiply(const float4x4& m, const float4& v);
  JTKQBVHDEF float4x4 matrix_matrix_multiply(const float4x4& left, const float4x4& right);
  JTKQBVHDEF float4x4 operator + (const float4x4& left, const float4x4& right);
  JTKQBVHDEF float4x4 operator - (const float4x4& left, const float4x4& right);
  JTKQBVHDEF float4x4 operator / (const float4x4& left, float value);
  JTKQBVHDEF float4x4 operator * (const float4x4& left, float value);
  JTKQBVHDEF float4x4 operator * (float value, const float4x4& right);


  struct hit
    {
    float u, v; // barycentric coordinates of hit
    float distance;
    int32_t found;
    };

  struct hit4
    {
    float4 u, v; // barycentric coordinates of hit
    float4 distance;
    bool4 found;
    };

  struct spherehit
    {
    float distance;
    int32_t found;
    };

  struct spherehit4
    {
    float4 distance;
    bool4 found;
    };

  struct ray
    {
    float4 orig;
    float4 dir;
    float t_near;
    float t_far;
    };

  struct woop_triangle
    {
    vec3<float4> v0;
    vec3<float4> v1;
    vec3<float4> v2;
    };

  struct woop_precompute
    {
    uint32_t kx, ky, kz;
    float4 Sx, Sy, Sz;
    };

  JTKQBVHDEF woop_precompute intersect_woop_precompute(const float4& r_dir);
  JTKQBVHDEF void intersect_woop(const woop_triangle& acc, const woop_precompute& pre, const vec3<float4>& r_orig, const float4& t_near, const float4& t_far, hit4& h);
  JTKQBVHDEF bool4 intersect(const float4* aabb, const vec3<float4>& r_origin, const float4& t_near, const float4& t_far, const int32_t* ray_dir_sign, const vec3<float4>& ray_inverse_dir);

  struct distance4
    {
    float4 u, v; // barycentric coordinates
    float4 distance_sqr;
    };

  JTKQBVHDEF void distance_sqr(const woop_triangle& acc, const vec3<float4>& point, distance4& dist);
  JTKQBVHDEF float4 distance_sqr(const float4* aabb, const vec3<float4>& point);

  JTKQBVHDEF float4x4 make_identity();
  JTKQBVHDEF vec3<float> transform(const float4x4& matrix, const vec3<float>& pt);
  JTKQBVHDEF vec3<float> transform_vector(const float4x4& matrix, const vec3<float>& vec);
  JTKQBVHDEF vec3<float> transform(const float4x4& matrix, const vec3<float>& pt, bool is_vector);
  JTKQBVHDEF float4 transform(const float4x4& matrix, const float4& pt);
  JTKQBVHDEF float4x4 make_transformation(const vec3<float>& i_origin, const vec3<float>& i_x_axis, const vec3<float>& i_y_axis, const vec3<float>& i_z_axis);
  JTKQBVHDEF float4x4 make_rotation(const vec3<float>& i_position, const vec3<float>& i_direction, float i_angle_radians);
  JTKQBVHDEF float4x4 make_scale3d(float scale_x, float scale_y, float scale_z);
  JTKQBVHDEF float4x4 make_translation(const vec3<float>& i_translation);
  JTKQBVHDEF vec3<float> get_translation(const float4x4& matrix);
  JTKQBVHDEF void set_x_axis(float4x4& matrix, const vec3<float>& x);
  JTKQBVHDEF void set_y_axis(float4x4& matrix, const vec3<float>& y);
  JTKQBVHDEF void set_z_axis(float4x4& matrix, const vec3<float>& z);
  JTKQBVHDEF void set_translation(float4x4& matrix, const vec3<float>& t);
  JTKQBVHDEF vec3<float> get_x_axis(const float4x4& matrix);
  JTKQBVHDEF vec3<float> get_y_axis(const float4x4& matrix);
  JTKQBVHDEF vec3<float> get_z_axis(const float4x4& matrix);
  JTKQBVHDEF float determinant(const float4x4& m);
  JTKQBVHDEF void solve_roll_pitch_yaw_transformation(float& rx, float& ry, float& rz, float& tx, float& ty, float& tz, const float4x4& m);
  JTKQBVHDEF float4x4 compute_from_roll_pitch_yaw_transformation(float rx, float ry, float rz, float tx, float ty, float tz);
  JTKQBVHDEF float4 roll_pitch_yaw_to_quaternion(float rx, float ry, float rz);
  JTKQBVHDEF float4x4 quaternion_to_rotation(const float4& quaternion);

  template <typename T>
  struct range
    {
    range();
    range(const T& begin, const T& end);
    T size() const;
    bool empty() const;
    T begin() const;
    T end() const;
    range intersect(const range& r) const;

    T _begin, _end;
    };

  template <class iterator, class predicate>
  iterator parallel_partition(iterator first, iterator last, predicate pred);


  struct qbvh_voxel
    {
    float4 bbox_min;
    float4 bbox_max;
    float4 centroid;

    void* operator new(size_t size); // for aligned allocation
    void operator delete(void* ptr); // for aligned allocation
    void* operator new[](size_t size); // for aligned allocation
    void operator delete[](void* ptr); // for aligned allocation
    };

  JTKQBVHDEF qbvh_voxel* build_triangle_qbvh_voxels(qbvh_voxel& total_bb, qbvh_voxel& centroid_bb, const vec3<float>* vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles);
  JTKQBVHDEF qbvh_voxel* build_sphere_qbvh_voxels(qbvh_voxel& total_bb, qbvh_voxel& centroid_bb, const vec3<float>* origins, const float* radii, uint32_t nr_of_spheres);

  struct qbvh_node
    {
    float4 bbox[2 * 3]; // 96 bytes  
    int4 child; // 16 bytes  
    int32_t free; // 4 bytes
    uint16_t nr_of_primitives[4]; // 8 bytes
    uint8_t axis0, axis1, axis2, pad; // 4 bytes 
    }; // 128 byte size
  JTKQBVHDEF void unite_four_aabbs(float4* out, float4* in, int k);
  JTKQBVHDEF void get_bbox(float4* bbox, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last, int k);
  JTKQBVHDEF float calculate_half_surface_area(const float4& left, const float4& right);
  JTKQBVHDEF uint8_t find_largest_dimension(const qbvh_voxel& bb);
  JTKQBVHDEF std::vector<uint32_t>::iterator partition(qbvh_voxel& bbox_left, qbvh_voxel& bbox_right, qbvh_voxel& centroid_left, qbvh_voxel& centroid_right, const uint32_t dim, const float split_pos, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);
  template <int K>
  void sah_optimized(std::vector<uint32_t>::iterator& mid, qbvh_voxel& bbox_left, qbvh_voxel& bbox_right, qbvh_voxel& centroid_left, qbvh_voxel& centroid_right, uint8_t& dim, const qbvh_voxel& bbox, const qbvh_voxel& centroid_bb, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);
  template <int K>
  void sah_parallel(std::vector<uint32_t>::iterator& mid, qbvh_voxel& bbox_left, qbvh_voxel& bbox_right, qbvh_voxel& centroid_left, qbvh_voxel& centroid_right, uint8_t& dim, const qbvh_voxel& bbox, const qbvh_voxel& centroid_bb, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

  class qbvh
    {
    public:

      struct properties
        {
        uint32_t leaf_size;
        int parallel_level;
        properties() : leaf_size(32), parallel_level(2) {}
        };

      qbvh(const std::vector<vec3<uint32_t>>& triangles, const vec3<float>* vertices);
      qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const std::vector<uint32_t>& ids, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr = properties());
      qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr = properties());
      hit find_closest_triangle(uint32_t& triangle_id, ray r, const vec3<uint32_t>* triangles, const vec3<float>* vertices) const;
      std::vector<hit> find_all_triangles(std::vector<uint32_t>& triangle_ids, ray r, const vec3<uint32_t>* triangles, const vec3<float>* vertices) const;
      hit find_closest_triangle(uint32_t& triangle_id, const vec3<float>& point, const vec3<uint32_t>* triangles, const vec3<float>* vertices) const;

    private:
      void _build(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb);
      typename int4::value_type construct_tree_single(aligned_vector<qbvh_node>& local_nodes, uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

      struct tree_stack
        {
        qbvh_voxel bb_child[4], centroid_child[4];
        uint32_t sizes[4];
        typename int4::value_type child[4];
        int level;
        std::vector<uint32_t>::iterator first;
        std::vector<uint32_t>::iterator last;
        std::vector<uint32_t>::iterator mid0, mid1, mid2;
        };

      typename int4::value_type construct_tree_prep(std::vector<tree_stack>& stack, uint32_t& sz, uint32_t& stack_index, uint32_t& node_index, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last, int level, bool prep);
      typename int4::value_type construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

    public:
      aligned_vector<qbvh_node> nodes;
      uint32_t root_id;
      std::vector<uint32_t>::iterator start;
      const properties props;
      std::vector<uint32_t> _ids;
    };

  JTKQBVHDEF void intersect_sphere(const vec3<float4>& sphere_origin, const float4& sphere_radius, const vec3<float4>& r_orig, const vec3<float4>& ray_dir, const float4& t_near, const float4& t_far, spherehit4& h);

  class sphere_qbvh
    {
    public:

      struct properties
        {
        uint32_t leaf_size;
        int parallel_level;
        properties() : leaf_size(32), parallel_level(2) {}
        };

      sphere_qbvh(const vec3<float>* origins, const float* radii, uint32_t nr_of_spheres);
      sphere_qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const std::vector<uint32_t>& ids, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr = properties());
      sphere_qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr = properties());
      spherehit find_closest_sphere(uint32_t& sphere_id, ray r, const vec3<float>* origins, const float* radii) const;
      std::vector<spherehit> find_all_spheres(std::vector<uint32_t>& sphere_ids, ray r, const vec3<float>* origins, const float* radii) const;

    private:
      void _build(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb);
      typename int4::value_type construct_tree_single(aligned_vector<qbvh_node>& local_nodes, uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

      struct tree_stack
        {
        qbvh_voxel bb_child[4], centroid_child[4];
        uint32_t sizes[4];
        typename int4::value_type child[4];
        int level;
        std::vector<uint32_t>::iterator first;
        std::vector<uint32_t>::iterator last;
        std::vector<uint32_t>::iterator mid0, mid1, mid2;
        };

      typename int4::value_type construct_tree_prep(std::vector<tree_stack>& stack, uint32_t& sz, uint32_t& stack_index, uint32_t& node_index, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last, int level, bool prep);
      typename int4::value_type construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

    public:
      aligned_vector<qbvh_node> nodes;
      uint32_t root_id;
      std::vector<uint32_t>::iterator start;
      const properties props;
      std::vector<uint32_t> _ids;
    };

  struct qbvh_two_level_node
    {
    float4 bbox[2 * 3]; // 96 bytes  
    int4 child; // 16 bytes  
    int32_t free1; // 4 bytes
    uint32_t free2; // 4 bytes
    uint32_t free3; // 4 bytes
    uint8_t axis0, axis1, axis2, pad; // 4 bytes  
    };

  JTKQBVHDEF void unite_four_aabbs(qbvh_voxel& out, const float4* in);
  class qbvh_two_level
    {
    enum
      {
      leaf_size = 1
      };

    public:
      qbvh_two_level(const qbvh** objects, uint32_t nr_of_objects);
      hit find_closest_triangle(uint32_t& triangle_id, uint32_t& object_id, ray r, const qbvh** objects, const vec3<uint32_t>** triangles, const vec3<float>** vertices) const;

    private:
      int32_t construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

    private:
      std::vector<uint32_t> object_ids;
      aligned_vector<qbvh_two_level_node> nodes;
      int32_t root_id;
      std::vector<uint32_t>::iterator start;
    };

  JTKQBVHDEF void transform(qbvh_voxel& voxel, const float4x4& matrix);
  class qbvh_two_level_with_transformations
    {
    enum
      {
      leaf_size = 1
      };

    public:
      qbvh_two_level_with_transformations(const qbvh** objects, const float4x4* transformations, uint32_t nr_of_objects);
      hit find_closest_triangle(uint32_t& triangle_id, uint32_t& object_id, ray r, const qbvh** objects, const float4x4* inverted_transformations, const vec3<uint32_t>** triangles, const vec3<float>** vertices) const;

    private:
      int32_t construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last);

    private:
      std::vector<uint32_t> object_ids;
      aligned_vector<qbvh_two_level_node> nodes;
      int32_t root_id;
      std::vector<uint32_t>::iterator start;
    };
    
    
  /////////////////////////////////////////////////////////////////////////
  // template class implementations
  /////////////////////////////////////////////////////////////////////////
    
  template<typename T, size_t alignment>
  typename aligned_allocator<T, alignment>::pointer aligned_allocator<T, alignment>::allocate(size_type n)
    {
    return (pointer)aligned_malloc(n * sizeof(value_type), alignment);
    }

  template<typename T, size_t alignment>
  void aligned_allocator<T, alignment>::deallocate(pointer p, size_type)
    {
    return aligned_free(p);
    }

  template<typename T, size_t alignment>
  void aligned_allocator<T, alignment>::construct(pointer p, const_reference val)
    {
    new (p) T(val);
    }

  template<typename T, size_t alignment>
  void aligned_allocator<T, alignment>::destroy(pointer p)
    {
    p->~T();
    }

  /////////////////////////////////////////////////////////////////////////
  // aligned vector
  /////////////////////////////////////////////////////////////////////////


  template<typename T, typename allocator>
  vector_t<T, allocator>::vector_t()
    : size_active(0), size_alloced(0), items(nullptr)
    {
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>::vector_t(size_t sz)
    : size_active(0), size_alloced(0), items(nullptr)
    {
    internal_resize_init(sz);
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>::vector_t(size_t sz, const allocator& alloc)
    : alloc(alloc), size_active(0), size_alloced(0), items(nullptr)
    {
    internal_resize_init(sz);
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>::vector_t(size_t sz, const T& val, const allocator& alloc)
    : alloc(alloc), size_active(0), size_alloced(0), items(nullptr)
    {
    internal_resize_fill(sz, val);
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>::~vector_t()
    {
    clear();
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>::vector_t(const vector_t& other)
    {
    size_active = other.size_active;
    size_alloced = other.size_alloced;
    items = alloc.allocate(size_alloced);
    for (size_t i = 0; i < size_active; i++)
      ::new (&items[i]) value_type(other.items[i]);
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>::vector_t(vector_t&& other)
    : alloc(std::move(other.alloc))
    {
    size_active = other.size_active; other.size_active = 0;
    size_alloced = other.size_alloced; other.size_alloced = 0;
    items = other.items; other.items = nullptr;
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>& vector_t<T, allocator>::operator=(const vector_t& other)
    {
    resize(other.size_active);
    for (size_t i = 0; i < size_active; i++)
      ::new (&items[i]) value_type(other.items[i]);
    return *this;
    }

  template<typename T, typename allocator>
  vector_t<T, allocator>& vector_t<T, allocator>::operator=(vector_t&& other)
    {
    vector_t temp(std::move(other));
    swap(temp);
    return *this;
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::swap(vector_t& other)
    {
    std::swap(alloc, other.alloc);
    std::swap(size_active, other.size_active);
    std::swap(size_alloced, other.size_alloced);
    std::swap(items, other.items);
    }

  template<typename T, typename allocator>
  typename vector_t<T, allocator>::iterator vector_t<T, allocator>::begin() { return items; };
  template<typename T, typename allocator>
  typename vector_t<T, allocator>::const_iterator vector_t<T, allocator>::begin() const { return items; };

  template<typename T, typename allocator>
  typename vector_t<T, allocator>::iterator vector_t<T, allocator>::end() { return items + size_active; };
  template<typename T, typename allocator>
  typename vector_t<T, allocator>::const_iterator vector_t<T, allocator>::end() const { return items + size_active; };

  template<typename T, typename allocator>
  bool vector_t<T, allocator>::empty() const { return size_active == 0; }
  template<typename T, typename allocator>
  size_t vector_t<T, allocator>::size() const { return size_active; }
  template<typename T, typename allocator>
  size_t vector_t<T, allocator>::capacity() const { return size_alloced; }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::resize(size_t new_size) {
    internal_resize(new_size, internal_grow_size(new_size));
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::reserve(size_t new_alloced)
    {
    /* do nothing if container already large enough */
    if (new_alloced <= size_alloced)
      return;

    /* resize exact otherwise */
    internal_resize(size_active, new_alloced);
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::shrink_to_fit() {
    internal_resize(size_active, size_active);
    }

  template<typename T, typename allocator>
   T& vector_t<T, allocator>::operator[](size_t i) { assert(i < size_active); return items[i]; }
  template<typename T, typename allocator>
  const T& vector_t<T, allocator>::operator[](size_t i) const { assert(i < size_active); return items[i]; }

  template<typename T, typename allocator>
        T& vector_t<T, allocator>::at(size_t i) { assert(i < size_active); return items[i]; }
  template<typename T, typename allocator>
  const T& vector_t<T, allocator>::at(size_t i) const { assert(i < size_active); return items[i]; }

  template<typename T, typename allocator>
  T& vector_t<T, allocator>::front() const { assert(size_active > 0); return items[0]; };
  template<typename T, typename allocator>
  T& vector_t<T, allocator>::back() const { assert(size_active > 0); return items[size_active - 1]; };

  template<typename T, typename allocator>
        T* vector_t<T, allocator>::data() { return items; };
  template<typename T, typename allocator>
  const T* vector_t<T, allocator>::data() const { return items; };

  template<typename T, typename allocator>
  void vector_t<T, allocator>::push_back(const T& nt)
    {
    const T v = nt; // need local copy as input reference could point to this vector
    internal_resize(size_active, internal_grow_size(size_active + 1));
    ::new (&items[size_active++]) T(v);
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::pop_back()
    {
    assert(!empty());
    size_active--;
    alloc.destroy(&items[size_active]);
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::clear()
    {
    /* destroy elements */
    for (size_t i = 0; i < size_active; i++)
      alloc.destroy(&items[i]);

    /* free memory */
    alloc.deallocate(items, size_alloced);
    items = nullptr;
    size_active = size_alloced = 0;
    }

  template<typename T, typename allocator>
  bool operator== (const vector_t<T, allocator>& a, const vector_t<T, allocator>& b)
    {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
      if (a[i] != b[i])
        return false;
    return true;
    }

  template<typename T, typename allocator>
  bool operator!= (const vector_t<T, allocator>& a, const vector_t<T, allocator>& b) {
    return !(a == b);
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::internal_resize_init(size_t new_active)
    {
    assert(size_active == 0);
    assert(size_alloced == 0);
    assert(items == nullptr);
    if (new_active == 0) return;
    items = alloc.allocate(new_active);
    for (size_t i = 0; i < new_active; i++) ::new (&items[i]) T();
    size_active = new_active;
    size_alloced = new_active;
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::internal_resize_fill(size_t new_active, const T& val)
    {
    assert(size_active == 0);
    assert(size_alloced == 0);
    assert(items == nullptr);
    if (new_active == 0) return;
    items = alloc.allocate(new_active);
    for (size_t i = 0; i < new_active; i++) ::new (&items[i]) T(val);
    size_active = new_active;
    size_alloced = new_active;
    }

  template<typename T, typename allocator>
  void vector_t<T, allocator>::internal_resize(size_t new_active, size_t new_alloced)
    {
    assert(new_active <= new_alloced);

    /* destroy elements */
    if (new_active < size_active)
      {
      for (size_t i = new_active; i < size_active; i++)
        alloc.destroy(&items[i]);
      size_active = new_active;
      }

    /* only reallocate if necessary */
    if (new_alloced == size_alloced) {
      for (size_t i = size_active; i < new_active; i++) ::new (&items[i]) T;
      size_active = new_active;
      return;
      }

    /* reallocate and copy items */
    T* old_items = items;
    items = alloc.allocate(new_alloced);
    for (size_t i = 0; i < size_active; i++) {
      ::new (&items[i]) T(std::move(old_items[i]));
      alloc.destroy(&old_items[i]);
      }
    for (size_t i = size_active; i < new_active; i++) {
      ::new (&items[i]) T;
      }
    alloc.deallocate(old_items, size_alloced);
    size_active = new_active;
    size_alloced = new_alloced;
    }

  template<typename T, typename allocator>
  size_t vector_t<T, allocator>::internal_grow_size(size_t new_alloced)
    {
    /* do nothing if container already large enough */
    if (new_alloced <= size_alloced)
      return size_alloced;

    /* resize to next power of 2 otherwise */
    size_t new_size_alloced = size_alloced;
    while (new_size_alloced < new_alloced) {
      new_size_alloced = 2 * new_size_alloced;
      if (new_size_alloced == 0) new_size_alloced = 1;
      }
    return new_size_alloced;
    }
    
  /////////////////////////////////////////////////////////////////////////
  // struct bool4
  /////////////////////////////////////////////////////////////////////////
    
  template <typename T>
  int32_t& bool4::operator [](T index) { return i[index]; }

  template <typename T>
  int32_t bool4::operator [](T index) const { return i[index]; }
  
  /////////////////////////////////////////////////////////////////////////
  // struct int4
  /////////////////////////////////////////////////////////////////////////
  
  template <typename T>
  int32_t& int4::operator [] (T ind)
    {
    return i[ind];
    }

  template <typename T>
  int32_t int4::operator [] (T ind) const
    {
    return i[ind];
    }
    
  /////////////////////////////////////////////////////////////////////////
  // struct float4
  /////////////////////////////////////////////////////////////////////////


  template <typename T>
  float& float4::operator [] (T i)
    {
    return f[i];
    }

  template <typename T>
  float float4::operator [] (T i) const
    {
    return f[i];
    }
    
  /////////////////////////////////////////////////////////////////////////
  // struct float4x4
  /////////////////////////////////////////////////////////////////////////

  // COLUMN MAJOR 4x4 MATRIX

  template <typename T>
  float& float4x4::operator [] (T i)
    {
    return f[i];
    }

  template <typename T>
  float float4x4::operator [] (T i) const
    {
    return f[i];
    }
    
  /////////////////////////////////////////////////////////////////////////
  // parallel partition
  /////////////////////////////////////////////////////////////////////////


  /*

From stack overflow:

I'd treat it as a degenerate case of parallel sample sort. (Parallel code for sample sort can be found here.)
Let N be the number of items. The degenerate sample sort will require O(N) temporary space, has
O(N) work, and O(P+ lg N) span (critical path). The last two values are important for analysis,
since speedup is limited to work/span.

I'm assuming the input is a random-access sequence. The steps are:

1. Allocate a temporary array big enough to hold a copy of the input sequence.

2. Divide the input into K blocks. K is a tuning parameter. For a system with P hardware threads,
K=max(4*P,L) might be good, where L is a constant for avoiding ridiculously small blocks.
The "4*P" allows some load balancing.

3. Move each block to its corresponding position in the temporary array and partition it using std::partition.
Blocks can be processed in parallel. Remember the offset of the "middle" for each block.
You might want to consider writing a custom routine that both moves (in the C++11 sense) and partitions a block.

4. Compute the offset to where each part of a block should go in the final result. The offsets for
the first part of each block can be done using an exclusive prefix sum over the offsets of the middles
from step 3. The offsets for the second part of each block can be computed similarly by using the
offset of each middle relative to the end of its block. The running sums in the latter case
become offsets from the end of the final output sequence. Unless you're dealing with more
than 100 hardware threads, I recommend using a serial exclusive scan.

5. Move the two parts of each block from the temporary array back to the appropriate places in the original sequence.
Copying each block can be done in parallel.

There is a way to embed the scan of step 4 into steps 3 and 5, so that the span can be
reduced to O(lg N), but I doubt it's worth the additional complexity.

If using tbb::parallel_for loops to parallelize steps 3 and 5, consider using
affinity_partitioner to help threads in step 5 pick up what they left in cache from step 3.

Note that partitioning requires only O(N) work for O(N) memory loads and stores.
Memory bandwidth could easily become the limiting resource for speedup.


JanM's notes: The extra array from step 1 kills all performance gain due to extra memory being used. Therefore
I'm following the same algorithm steps, but do everything in place.


*/

  template <typename T>
  range<T>::range() {}

  template <typename T>
  range<T>::range(const T& begin, const T& end) : _begin(begin), _end(end) {}

  template <typename T>
  T range<T>::size() const { return _end - _begin; }

  template <typename T>
  bool range<T>::empty() const { return _end <= _begin; }
  template <typename T>
  T range<T>::begin() const { return _begin; }
  template <typename T>
  T range<T>::end() const { return _end; }
  template <typename T>
  range<T> range<T>::intersect(const range<T>& r) const
    {
    return range(std::max(_begin, r._begin), std::min(_end, r._end));
    }

  namespace parallel_partition_details
    {
    inline const range<uint64_t>* findStartRange(uint64_t& index, const range<uint64_t>* const r)
      {
      uint64_t i = 0;
      while (index >= (uint64_t)r[i].size())
        {
        index -= (uint64_t)r[i].size();
        i++;
        }
      return &r[i];
      }
    }

  template <class iterator, class predicate>
  iterator parallel_partition(iterator first, iterator last, predicate pred)
    {
    const auto N = last - first;

    if (N < 16384)
      return std::partition(first, last, pred);

    const unsigned int P = hardware_concurrency();
    const unsigned int K = std::min<unsigned int>(std::max<unsigned int>(P << 2, 16), 64);

    if (P == 1)
      return std::partition(first, last, pred);

#ifdef _WIN32
    _declspec(align(64))
#endif
      uint64_t offset_from_first_in_block[64]
#ifndef _WIN32 // linux alignment in gcc
      __attribute__((aligned(64)))
#endif
      ;

#ifdef _WIN32
    _declspec(align(64))
#endif
      uint64_t offset_from_last_in_block[64]
#ifndef _WIN32 // linux alignment in gcc
      __attribute__((aligned(64)))
#endif
      ;

    parallel_for((unsigned int)0, K, [&](unsigned int k)
      {
      const auto s = (uint64_t)k * N / (uint64_t)K;
      const auto e = (uint64_t)(k + 1) * N / (uint64_t)K;
      const auto sz = e - s;
      auto block_first = first + s;
      auto block_last = block_first + sz;
      auto mid = std::partition(block_first, block_last, pred);
      offset_from_first_in_block[k] = mid - block_first;
      offset_from_last_in_block[k] = sz - offset_from_first_in_block[k];
      });
#ifdef _WIN32
    _declspec(align(64))
#endif
      uint64_t offset_from_first[65]
#ifndef _WIN32 // linux alignment in gcc
      __attribute__((aligned(64)))
#endif
      ;
#ifdef _WIN32
    _declspec(align(64))
#endif
      uint64_t offset_from_last[65]
#ifndef _WIN32 // linux alignment in gcc
      __attribute__((aligned(64)))
#endif
      ;
    offset_from_first[0] = 0;
    offset_from_last[0] = 0;
    for (unsigned int i = 0; i < K; ++i)
      {
      offset_from_first[i + 1] = offset_from_first[i] + offset_from_first_in_block[i];
      offset_from_last[i + 1] = offset_from_last[i] + offset_from_last_in_block[i];
      }

    const auto mid = first + offset_from_first[K];

    const range<uint64_t> global_left(0, offset_from_first[K]);
    const range<uint64_t> global_right(offset_from_first[K], N);

    uint64_t numMisplacedRangesLeft = 0;
    uint64_t numMisplacedRangesRight = 0;
    uint64_t numMisplacedItemsLeft = 0;
    uint64_t numMisplacedItemsRight = 0;

#ifdef _WIN32
    _declspec(align(64))
#endif
      range<uint64_t> leftMisplacedRanges[64]
#ifndef _WIN32 // linux alignment in gcc
      __attribute__((aligned(64)))
#endif
      ;
#ifdef _WIN32
    _declspec(align(64))
#endif
      range<uint64_t> rightMisplacedRanges[64]
#ifndef _WIN32 // linux alignment in gcc
      __attribute__((aligned(64)))
#endif
      ;

    for (unsigned int k = 0; k < K; ++k)
      {
      const auto s = (uint64_t)k * N / (uint64_t)K;
      const auto e = (uint64_t)(k + 1) * N / (uint64_t)K;
      const range<uint64_t> left_range(s, s + offset_from_first_in_block[k]);
      const range<uint64_t> right_range(s + offset_from_first_in_block[k], e);
      const range<uint64_t> left_misplaced = global_left.intersect(right_range);
      const range<uint64_t> right_misplaced = global_right.intersect(left_range);
      if (!left_misplaced.empty())
        {
        numMisplacedItemsLeft += left_misplaced.size();
        leftMisplacedRanges[numMisplacedRangesLeft++] = left_misplaced;
        }

      if (!right_misplaced.empty())
        {
        numMisplacedItemsRight += right_misplaced.size();
        rightMisplacedRanges[numMisplacedRangesRight++] = right_misplaced;
        }
      }
    assert(numMisplacedItemsLeft == numMisplacedItemsRight);

    if (numMisplacedItemsLeft == 0)
      return mid;

    parallel_for((unsigned int)0, K, [&](unsigned int k)
      {
      const auto s = (uint64_t)k * numMisplacedItemsLeft / (uint64_t)K;
      const auto e = (uint64_t)(k + 1) * numMisplacedItemsLeft / (uint64_t)K;
      uint64_t leftLocalIndex = s;
      uint64_t rightLocalIndex = s;
      const range<uint64_t>* l_range = parallel_partition_details::findStartRange(leftLocalIndex, leftMisplacedRanges);
      const range<uint64_t>* r_range = parallel_partition_details::findStartRange(rightLocalIndex, rightMisplacedRanges);
      uint64_t l_left = l_range->size() - leftLocalIndex;
      uint64_t r_left = r_range->size() - rightLocalIndex;
      iterator l = first + l_range->begin() + leftLocalIndex;
      iterator r = first + r_range->begin() + rightLocalIndex;
      uint64_t size = e - s;
      uint64_t items = std::min<uint64_t>(size, std::min<uint64_t>(l_left, r_left));

      while (size)
        {
        if (l_left == 0)
          {
          ++l_range;
          l_left = l_range->size();
          l = first + l_range->begin();
          items = std::min<uint64_t>(size, std::min<uint64_t>(l_left, r_left));
          }

        if (r_left == 0)
          {
          ++r_range;
          r_left = r_range->size();
          r = first + r_range->begin();
          items = std::min<uint64_t>(size, std::min<uint64_t>(l_left, r_left));
          }

        size -= items;
        l_left -= items;
        r_left -= items;

        while (items)
          {
          --items;
          const typename iterator::value_type tmp = *l;
          *l++ = *r;
          *r++ = tmp;
          }
        }
      });

    return mid;
    }
    
  static const uint32_t order[8][4] =
    { { 3,2,1,0 },
    { 3,2,0,1 },
    { 2,3,1,0 },
    { 2,3,0,1 },
    { 1,0,3,2 },
    { 1,0,2,3 },
    { 0,1,3,2 },
    { 0,1,2,3 }
    };


  static const uint32_t sign_bit = 0x80000000;
  static const int4 sign_mask = int4(0x80000000, 0x80000000, 0x80000000, 0x80000000);
  static const int4 bit_mask = int4(0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF);
  
  static const __m128 lookup_mask[16] =
    {
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)), // 0000b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0xffffffff)), // 0001b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0x00000000, 0xffffffff, 0x00000000)), // 0010b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0x00000000, 0xffffffff, 0xffffffff)), // 0011b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0xffffffff, 0x00000000, 0x00000000)), // 0100b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0xffffffff, 0x00000000, 0xffffffff)), // 0101b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0xffffffff, 0xffffffff, 0x00000000)), // 0110b
    _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0xffffffff, 0xffffffff, 0xffffffff)), // 0111b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0x00000000, 0x00000000, 0x00000000)), // 1000b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0x00000000, 0x00000000, 0xffffffff)), // 1001b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0x00000000, 0xffffffff, 0x00000000)), // 1010b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0x00000000, 0xffffffff, 0xffffffff)), // 1011b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0xffffffff, 0x00000000, 0x00000000)), // 1100b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0xffffffff, 0x00000000, 0xffffffff)), // 1101b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0x00000000)), // 1110b
    _mm_castsi128_ps(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff))  // 1111b
    };

  static const float epsilon = std::numeric_limits<float>::epsilon()*2.f;
  static const __m128 three_128 = _mm_set1_ps(3.0f);
  static const __m128 half_128 = _mm_set1_ps(0.5f);
  static const __m128 one_128 = _mm_set1_ps(1.f);
  static const __m128 zero_128 = _mm_set1_ps(0.f);
  static const __m128i not_zero = _mm_set_epi32(~0, ~0, ~0, ~0);
  static const __m128 epsilon_128 = _mm_set1_ps(epsilon);
  static const __m128 minus_epsilon_128 = _mm_set1_ps(-epsilon);
  static const __m128 one_plus_epsilon_128 = _mm_set1_ps(1.f + epsilon);
  static const uint32_t modulo[] = { 0,1,2,0,1 };
    
  JTKQBVHINLINE bool4::bool4() {}
  JTKQBVHINLINE bool4::bool4(const __m128 in) : m128(in) {}
  JTKQBVHINLINE bool4::bool4(const __m128i in) : m128(_mm_castsi128_ps(in)) {}
  JTKQBVHINLINE bool4::bool4(bool b) : m128(lookup_mask[b ? 15 : 0]) {}
  JTKQBVHINLINE bool4::bool4(bool b0, bool b1, bool b2, bool b3) : m128(lookup_mask[(uint8_t(b3) << 3) | (uint8_t(b2) << 2) | (uint8_t(b1) << 1) | (uint8_t(b0))]) {}

  JTKQBVHINLINE int4::int4() {}
  JTKQBVHINLINE int4::int4(const __m128i& in) : m128i(in) {}
  JTKQBVHINLINE int4::int4(int32_t in) : m128i(_mm_set1_epi32(in)) {}
  JTKQBVHINLINE int4::int4(int32_t i0, int32_t i1, int32_t i2, int32_t i3) : m128i(_mm_set_epi32(i3, i2, i1, i0)) {}
  JTKQBVHINLINE int4::int4(const bool4& in) : m128i(_mm_castps_si128(in.m128)) {}
  
  JTKQBVHINLINE float4::float4() {}
  JTKQBVHINLINE float4::float4(const __m128 in) : m128(in) {}
  JTKQBVHINLINE float4::float4(float f) : m128(_mm_set1_ps(f)) {}
  JTKQBVHINLINE float4::float4(float _x, float _y, float _z) : m128(_mm_set_ps(1.f, _z, _y, _x)) {}
  JTKQBVHINLINE float4::float4(float _x, float _y, float _z, float _w) : m128(_mm_set_ps(_w, _z, _y, _x)) {}
  
  JTKQBVHINLINE float4x4::float4x4() {}
  JTKQBVHINLINE float4x4::float4x4(const float4& col0, const float4& col1, const float4& col2, const float4& col3) : col{ col0, col1, col2, col3 } {}
  JTKQBVHINLINE float4x4::float4x4(float* m)
    {
    for (int i = 0; i < 16; ++i)
      f[i] = m[i];
    }
    
  JTKQBVHINLINE void* qbvh_voxel::operator new(size_t size) { return aligned_malloc(size, 16); }
  JTKQBVHINLINE void qbvh_voxel::operator delete(void* ptr) { aligned_free(ptr); }
  JTKQBVHINLINE void* qbvh_voxel::operator new[](size_t size) { return aligned_malloc(size, 16); }
  JTKQBVHINLINE void qbvh_voxel::operator delete[](void* ptr) { aligned_free(ptr); }


  template <int K>
  void sah_optimized(std::vector<uint32_t>::iterator& mid, qbvh_voxel& bbox_left, qbvh_voxel& bbox_right, qbvh_voxel& centroid_left, qbvh_voxel& centroid_right, uint8_t& dim, const qbvh_voxel& bbox, const qbvh_voxel& centroid_bb, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    const uint32_t nr_of_triangles = (uint32_t)std::distance(first, last);
    const float eps = std::numeric_limits<float>::epsilon() * 1024;
    const float one_minus_eps = (float)1 - eps;
    const float one_over_K = (float)1 / (float)K;
    std::array<uint32_t, 2 * K> bin;
    memset(&bin[0], 0, sizeof(uint32_t) * 2 * K);

    dim = find_largest_dimension(centroid_bb);

    const float s = bbox.bbox_max[dim] - bbox.bbox_min[dim];
    const float k1 = (float)K * one_minus_eps / s;
    const float k0 = bbox.bbox_min[dim];


    for (auto it = first; it != last; ++it)
      {
      const uint32_t bmin = (uint32_t)std::floor(k1*(voxels[*it].bbox_min[dim] - k0));
      const uint32_t bmax = (uint32_t)std::floor(k1*(voxels[*it].bbox_max[dim] - k0));
      ++bin[bmin];
      ++bin[K + bmax];
      }

    const float total_surface_area = calculate_half_surface_area(bbox.bbox_min, bbox.bbox_max);

    float inv_total_surface_area = 0;

    if (total_surface_area > eps)
      inv_total_surface_area = (float)1 / total_surface_area;

    float split_pos;
    float minimal_cost;


    const float bsize = s;
    const float bstep = bsize * one_over_K;

    split_pos = k0 + (float)0.5*bstep;
    minimal_cost = std::numeric_limits<float>::max();

    uint32_t left = 0;
    uint32_t right = nr_of_triangles;
    float4 bminleft = bbox.bbox_min;
    float4 bminright = bbox.bbox_min;
    float4 bmaxleft = bbox.bbox_max;
    float4 bmaxright = bbox.bbox_max;

    for (int i = 0; i < K - 1; ++i)
      {
      left += bin[i];
      right -= bin[K + i];

      assert(left <= nr_of_triangles);
      assert(right <= nr_of_triangles);

      const float pos = k0 + (i + 0.5f)*bstep;
      bmaxleft[dim] = pos;
      bminright[dim] = pos;

      const float saleft = calculate_half_surface_area(bminleft, bmaxleft);
      const float saright = calculate_half_surface_area(bminright, bmaxright);

      const float twice_Taabb = (float)0.4;
      const float Ttri = (float)0.8;
      //float cost = 2.0*Taabb + (saleft * inv_total_surface_area) * left * Ttri + (saright * inv_total_surface_area) * right * Ttri;
      const float cost = twice_Taabb + (saleft * left + saright * right) * inv_total_surface_area * Ttri;

      if (cost < minimal_cost)
        {
        minimal_cost = cost;
        split_pos = pos;
        }
      }

    bbox_left.bbox_min = std::numeric_limits<float>::max();
    bbox_left.bbox_max = -std::numeric_limits<float>::max();

    bbox_right.bbox_min = bbox_left.bbox_min;
    bbox_right.bbox_max = bbox_left.bbox_max;

    centroid_left.bbox_min = bbox_left.bbox_min;
    centroid_left.bbox_max = bbox_left.bbox_max;

    centroid_right.bbox_min = bbox_left.bbox_min;
    centroid_right.bbox_max = bbox_left.bbox_max;

    mid = partition(bbox_left, bbox_right, centroid_left, centroid_right, dim, split_pos, voxels, first, last);
    if ((mid == first) || (mid == last))
      {
      bbox_left.bbox_min = std::numeric_limits<float>::max();
      bbox_left.bbox_max = -std::numeric_limits<float>::max();

      bbox_right.bbox_min = bbox_left.bbox_min;
      bbox_right.bbox_max = bbox_left.bbox_max;

      centroid_left.bbox_min = bbox_left.bbox_min;
      centroid_left.bbox_max = bbox_left.bbox_max;

      centroid_right.bbox_min = bbox_left.bbox_min;
      centroid_right.bbox_max = bbox_left.bbox_max;

      mid = first + (last - first) / 2;

      std::nth_element(first, mid, last, [&](uint32_t id1, uint32_t id2)
        {
        return voxels[id1].centroid[dim] < voxels[id2].centroid[dim];
        });

      for (auto it = first; it != mid; ++it)
        {
        bbox_left.bbox_min = min(bbox_left.bbox_min, voxels[*it].bbox_min);
        bbox_left.bbox_max = max(bbox_left.bbox_max, voxels[*it].bbox_max);
        centroid_left.bbox_min = min(centroid_left.bbox_min, voxels[*it].centroid);
        centroid_left.bbox_max = max(centroid_left.bbox_max, voxels[*it].centroid);
        }
      for (auto it = mid; it != last; ++it)
        {
        bbox_right.bbox_min = min(bbox_right.bbox_min, voxels[*it].bbox_min);
        bbox_right.bbox_max = max(bbox_right.bbox_max, voxels[*it].bbox_max);
        centroid_right.bbox_min = min(centroid_right.bbox_min, voxels[*it].centroid);
        centroid_right.bbox_max = max(centroid_right.bbox_max, voxels[*it].centroid);
        }
      }
    }


  template <int K>
  void sah_parallel(std::vector<uint32_t>::iterator& mid, qbvh_voxel& bbox_left, qbvh_voxel& bbox_right, qbvh_voxel& centroid_left, qbvh_voxel& centroid_right, uint8_t& dim, const qbvh_voxel& bbox, const qbvh_voxel& centroid_bb, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    const uint32_t nr_of_triangles = (uint32_t)std::distance(first, last);
    const float eps = std::numeric_limits<float>::epsilon() * 1024;
    const float one_minus_eps = (float)1 - eps;
    const float one_over_K = (float)1 / (float)K;
    std::array<uint32_t, 2 * K> bin;
    memset(&bin[0], 0, sizeof(uint32_t) * 2 * K);

    dim = find_largest_dimension(centroid_bb);

    const float sw = bbox.bbox_max[dim] - bbox.bbox_min[dim];
    const float k1 = (float)K * one_minus_eps / sw;
    const float k0 = bbox.bbox_min[dim];

    static const unsigned int number_of_threads = hardware_concurrency();
    static const unsigned int number_of_blocks = number_of_threads * 4;

    std::vector<std::array<uint32_t, 2 * K>> local_bins(number_of_blocks, bin);

    parallel_for((unsigned int)0, number_of_blocks, [&](unsigned int t)
      {
      const uint32_t s = (uint32_t)((uint64_t)t * (uint64_t)(nr_of_triangles) / (uint64_t)number_of_blocks);
      const uint32_t e = (uint32_t)((uint64_t)(t + 1) * (uint64_t)(nr_of_triangles) / (uint64_t)number_of_blocks);
      auto it = first + s;
      const auto it_end = first + e;
      for (; it != it_end; ++it)
        {
        const uint32_t bmin = (uint32_t)std::floor(k1*(voxels[*it].bbox_min[dim] - k0));
        const uint32_t bmax = (uint32_t)std::floor(k1*(voxels[*it].bbox_max[dim] - k0));
        ++local_bins[t][bmin];
        ++local_bins[t][K + bmax];
        }
      });

    parallel_for(uint32_t(0), uint32_t(2 * K), [&](uint32_t i)
      {
      for (unsigned int t = 0; t < number_of_blocks; ++t)
        bin[i] += local_bins[t][i];
      });

    const float total_surface_area = calculate_half_surface_area(bbox.bbox_min, bbox.bbox_max);

    float inv_total_surface_area = 0;

    if (total_surface_area > eps)
      inv_total_surface_area = (float)1 / total_surface_area;

    float split_pos;
    float minimal_cost;


    const float bsize = sw;
    const float bstep = bsize * one_over_K;

    split_pos = k0 + (float)0.5*bstep;
    minimal_cost = std::numeric_limits<float>::max();

    uint32_t left = 0;
    uint32_t right = nr_of_triangles;
    float4 bminleft = bbox.bbox_min;
    float4 bminright = bbox.bbox_min;
    float4 bmaxleft = bbox.bbox_max;
    float4 bmaxright = bbox.bbox_max;

    for (int i = 0; i < K - 1; ++i)
      {
      left += bin[i];
      right -= bin[K + i];

      assert(left <= nr_of_triangles);
      assert(right <= nr_of_triangles);

      const float pos = k0 + (i + (float)0.5)*bstep;
      bmaxleft[dim] = pos;
      bminright[dim] = pos;

      const float saleft = calculate_half_surface_area(bminleft, bmaxleft);
      const float saright = calculate_half_surface_area(bminright, bmaxright);

      const float twice_Taabb = (float)0.4;
      const float Ttri = (float)0.8;
      const float cost = twice_Taabb + (saleft * left + saright * right) * inv_total_surface_area * Ttri;

      if (cost < minimal_cost)
        {
        minimal_cost = cost;
        split_pos = pos;
        }
      }

    mid = parallel_partition(first, last, [&](uint32_t id)
      {
      return voxels[id].centroid[dim] < split_pos;
      });
    if ((mid == first) || (mid == last))
      {
      mid = first + (last - first) / 2;
      std::nth_element(first, mid, last, [&](uint32_t id1, uint32_t id2)
        {
        return voxels[id1].centroid[dim] < voxels[id2].centroid[dim];
        });
      }

    bbox_left.bbox_min = std::numeric_limits<float>::max();
    bbox_left.bbox_max = -std::numeric_limits<float>::max();

    bbox_right.bbox_min = bbox_left.bbox_min;
    bbox_right.bbox_max = bbox_left.bbox_max;

    centroid_left.bbox_min = bbox_left.bbox_min;
    centroid_left.bbox_max = bbox_left.bbox_max;

    centroid_right.bbox_min = bbox_left.bbox_min;
    centroid_right.bbox_max = bbox_left.bbox_max;

    aligned_vector<qbvh_voxel> bbox_left_blocks(number_of_blocks, bbox_left);
    aligned_vector<qbvh_voxel> bbox_right_blocks(number_of_blocks, bbox_left);
    aligned_vector<qbvh_voxel> centroid_left_blocks(number_of_blocks, bbox_left);
    aligned_vector<qbvh_voxel> centroid_right_blocks(number_of_blocks, bbox_left);

    const uint32_t nr_of_left_items = (uint32_t)std::distance(first, mid);
    const uint32_t nr_of_right_items = (uint32_t)std::distance(mid, last);

    parallel_for((unsigned int)0, number_of_blocks, [&](unsigned int t)
      {
      const uint64_t s1 = (uint64_t)t * (uint64_t)(nr_of_left_items) / (uint64_t)number_of_blocks;
      const uint64_t e1 = (uint64_t)(t + 1) *(uint64_t)(nr_of_left_items) / (uint64_t)number_of_blocks;
      auto it1 = first + s1;
      const auto it1_end = first + e1;
      for (; it1 != it1_end; ++it1)
        {
        bbox_left_blocks[t].bbox_min = min(bbox_left_blocks[t].bbox_min, voxels[*it1].bbox_min);
        bbox_left_blocks[t].bbox_max = max(bbox_left_blocks[t].bbox_max, voxels[*it1].bbox_max);
        centroid_left_blocks[t].bbox_min = min(centroid_left_blocks[t].bbox_min, voxels[*it1].centroid);
        centroid_left_blocks[t].bbox_max = max(centroid_left_blocks[t].bbox_max, voxels[*it1].centroid);
        }

      const uint32_t s2 = (uint64_t)t * (uint64_t)(nr_of_right_items) / (uint64_t)number_of_blocks;
      const uint32_t e2 = (uint64_t)(t + 1) * (uint64_t)(nr_of_right_items) / (uint64_t)number_of_blocks;
      auto it2 = mid + s2;
      const auto it2_end = mid + e2;
      for (; it2 != it2_end; ++it2)
        {
        bbox_right_blocks[t].bbox_min = min(bbox_right_blocks[t].bbox_min, voxels[*it2].bbox_min);
        bbox_right_blocks[t].bbox_max = max(bbox_right_blocks[t].bbox_max, voxels[*it2].bbox_max);
        centroid_right_blocks[t].bbox_min = min(centroid_right_blocks[t].bbox_min, voxels[*it2].centroid);
        centroid_right_blocks[t].bbox_max = max(centroid_right_blocks[t].bbox_max, voxels[*it2].centroid);
        }
      });

    for (unsigned int t = 0; t < number_of_blocks; ++t)
      {
      bbox_left.bbox_min = min(bbox_left.bbox_min, bbox_left_blocks[t].bbox_min);
      bbox_left.bbox_max = max(bbox_left.bbox_max, bbox_left_blocks[t].bbox_max);
      centroid_left.bbox_min = min(centroid_left.bbox_min, centroid_left_blocks[t].bbox_min);
      centroid_left.bbox_max = max(centroid_left.bbox_max, centroid_left_blocks[t].bbox_max);
      bbox_right.bbox_min = min(bbox_right.bbox_min, bbox_right_blocks[t].bbox_min);
      bbox_right.bbox_max = max(bbox_right.bbox_max, bbox_right_blocks[t].bbox_max);
      centroid_right.bbox_min = min(centroid_right.bbox_min, centroid_right_blocks[t].bbox_min);
      centroid_right.bbox_max = max(centroid_right.bbox_max, centroid_right_blocks[t].bbox_max);
      }
    }

  JTKQBVHINLINE qbvh::qbvh(const std::vector<vec3<uint32_t>>& triangles, const vec3<float>* vertices)
    {
    qbvh_voxel total_bb, centroid_bb;
    const uint32_t nr_of_triangles = (uint32_t)triangles.size();
    auto qbvh_voxels = build_triangle_qbvh_voxels(total_bb, centroid_bb, vertices, triangles.data(), nr_of_triangles);
    _build(qbvh_voxels, nr_of_triangles, total_bb, centroid_bb);
    delete[] qbvh_voxels;
    }

  JTKQBVHINLINE qbvh::qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const std::vector<uint32_t>& ids, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr) : props(pr)
    {
    assert(props.leaf_size > 1);
    _ids = ids;
    _build(voxels, nr_of_items, total_bb, centroid_bb);
    }

  JTKQBVHINLINE qbvh::qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr) : props(pr)
    {
    assert(props.leaf_size > 1);
    _build(voxels, nr_of_items, total_bb, centroid_bb);
    }

  JTKQBVHINLINE hit qbvh::find_closest_triangle(uint32_t& triangle_id, ray r, const vec3<uint32_t>* triangles, const vec3<float>* vertices) const
    {
    hit h;
    h.found = 0;
    h.distance = std::numeric_limits<float>::max();

    auto woop_precomputation = intersect_woop_precompute(r.dir);

    vec3<float4> ray_origin(r.orig[0], r.orig[1], r.orig[2]);
    vec3<float4> ray_dir(r.dir[0], r.dir[1], r.dir[2]);

    int32_t ray_dir_sign[3];
    float4 ray_inverse_dir_tmp = reciprocal(r.dir);
    if (fabs(ray_inverse_dir_tmp[0]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[0] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[1]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[1] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[2]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[2] = std::numeric_limits<float>::max();
    vec3<float4> ray_inverse_dir(ray_inverse_dir_tmp[0], ray_inverse_dir_tmp[1], ray_inverse_dir_tmp[2]);
    float4 t_near(r.t_near);
    float4 t_far(r.t_far);
    ray_dir_sign[0] = r.dir[0] < 0 ? 1 : 0;
    ray_dir_sign[1] = r.dir[1] < 0 ? 1 : 0;
    ray_dir_sign[2] = r.dir[2] < 0 ? 1 : 0;

    std::vector<int32_t> stack;
    stack.reserve(512);
    stack.push_back((int32_t)root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();

      const uint32_t order_pos = (ray_dir_sign[node.axis0] << 2) | (ray_dir_sign[node.axis1] << 1) | (ray_dir_sign[node.axis2]);
      const uint32_t* local_order = order[order_pos];

      auto isct = intersect(node.bbox, ray_origin, t_near, t_far, ray_dir_sign, ray_inverse_dir);

      if (any(isct))
        {
        const int4 is_leaf = node.child & sign_mask;
        const int contains_at_least_one_leaf = any(is_leaf);
        if (contains_at_least_one_leaf)
          {
          const int4 index = node.child & bit_mask;
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              {
              if (is_leaf[o])
                {
                const uint32_t d = ((uint32_t)node.nr_of_primitives[o] + 3) >> 2;
                if (d) // not empty
                  {
                  uint64_t ind = index[o];
                  for (uint16_t k = 0; k < d; ++k)
                    {
                    woop_triangle acc;
                    const vec3<uint32_t>* tria0 = triangles + _ids[ind];
                    const uint32_t v00 = (*tria0)[0];
                    const uint32_t v01 = (*tria0)[1];
                    const uint32_t v02 = (*tria0)[2];
                    const vec3<uint32_t>* tria1 = triangles + _ids[++ind];
                    const uint32_t v10 = (*tria1)[0];
                    const uint32_t v11 = (*tria1)[1];
                    const uint32_t v12 = (*tria1)[2];
                    const vec3<uint32_t>* tria2 = triangles + _ids[++ind];
                    const uint32_t v20 = (*tria2)[0];
                    const uint32_t v21 = (*tria2)[1];
                    const uint32_t v22 = (*tria2)[2];
                    const vec3<uint32_t>* tria3 = triangles + _ids[++ind];
                    const uint32_t v30 = (*tria3)[0];
                    const uint32_t v31 = (*tria3)[1];
                    const uint32_t v32 = (*tria3)[2];
                    ++ind;

                    _mm_prefetch((char *)(vertices + v00), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v10), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v20), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v30), _MM_HINT_T0);

                    _mm_prefetch((char *)(vertices + v01), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v11), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v21), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v31), _MM_HINT_T0);

                    _mm_prefetch((char *)(vertices + v02), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v12), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v22), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v32), _MM_HINT_T0);

                    acc.v0[0] = float4(vertices[v00][0], vertices[v10][0], vertices[v20][0], vertices[v30][0]);
                    acc.v0[1] = float4(vertices[v00][1], vertices[v10][1], vertices[v20][1], vertices[v30][1]);
                    acc.v0[2] = float4(vertices[v00][2], vertices[v10][2], vertices[v20][2], vertices[v30][2]);


                    acc.v1[0] = float4(vertices[v01][0], vertices[v11][0], vertices[v21][0], vertices[v31][0]);
                    acc.v1[1] = float4(vertices[v01][1], vertices[v11][1], vertices[v21][1], vertices[v31][1]);
                    acc.v1[2] = float4(vertices[v01][2], vertices[v11][2], vertices[v21][2], vertices[v31][2]);


                    acc.v2[0] = float4(vertices[v02][0], vertices[v12][0], vertices[v22][0], vertices[v32][0]);
                    acc.v2[1] = float4(vertices[v02][1], vertices[v12][1], vertices[v22][1], vertices[v32][1]);
                    acc.v2[2] = float4(vertices[v02][2], vertices[v12][2], vertices[v22][2], vertices[v32][2]);

                    hit4 local_h;
                    intersect_woop(acc, woop_precomputation, ray_origin, t_near, t_far, local_h);
                    //intersect_moeller_trumbore(acc, ray_origin, ray_dir, t_near, t_far, local_h);
                    const bool4 mask = (abs(local_h.distance) < float4(std::abs(h.distance))) & local_h.found;
                    if (any(mask))
                      {
                      for (int m = 0; m < 4; ++m)
                        {
                        if (mask[m])
                          {
                          if (std::abs(local_h.distance[m]) < std::abs(h.distance))
                            {
                            h.found = -1;
                            h.distance = local_h.distance[m];
                            triangle_id = _ids[index[o] + k * 4 + m];
                            h.u = local_h.u[m];
                            h.v = local_h.v[m];
                            (h.distance > 0) ? t_far = float4(h.distance) : t_near = float4(h.distance);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              else
                stack.push_back((int32_t)node.child[o]);
              }
            }
          }
        else
          {
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              {
              stack.push_back((int32_t)node.child[o]);
              }
            }
          }
        }
      }
    return h;
    }

  JTKQBVHINLINE std::vector<hit> qbvh::find_all_triangles(std::vector<uint32_t>& triangle_ids, ray r, const vec3<uint32_t>* triangles, const vec3<float>* vertices) const
    {
    std::vector<hit> hits;

    auto woop_precomputation = intersect_woop_precompute(r.dir);

    vec3<float4> ray_origin(r.orig[0], r.orig[1], r.orig[2]);
    vec3<float4> ray_dir(r.dir[0], r.dir[1], r.dir[2]);

    int32_t ray_dir_sign[3];
    float4 ray_inverse_dir_tmp = reciprocal(r.dir);
    if (fabs(ray_inverse_dir_tmp[0]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[0] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[1]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[1] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[2]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[2] = std::numeric_limits<float>::max();
    vec3<float4> ray_inverse_dir(ray_inverse_dir_tmp[0], ray_inverse_dir_tmp[1], ray_inverse_dir_tmp[2]);
    float4 t_near(r.t_near);
    float4 t_far(r.t_far);
    ray_dir_sign[0] = r.dir[0] < 0 ? 1 : 0;
    ray_dir_sign[1] = r.dir[1] < 0 ? 1 : 0;
    ray_dir_sign[2] = r.dir[2] < 0 ? 1 : 0;

    std::vector<int32_t> stack;
    stack.reserve(512);
    stack.push_back((int32_t)root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();

      auto isct = intersect(node.bbox, ray_origin, t_near, t_far, ray_dir_sign, ray_inverse_dir);

      if (any(isct))
        {
        const int4 is_leaf = node.child & sign_mask;
        const int contains_at_least_one_leaf = any(is_leaf);
        if (contains_at_least_one_leaf)
          {
          const int4 index = node.child & bit_mask;
          for (int i = 0; i < 4; ++i)
            {
            if (isct[i])
              {
              if (is_leaf[i])
                {
                const uint32_t d = ((uint32_t)node.nr_of_primitives[i] + 3) >> 2;
                if (d) // not empty
                  {
                  uint64_t ind = index[i];
                  for (uint16_t k = 0; k < d; ++k)
                    {
                    woop_triangle acc;
                    const vec3<uint32_t>* tria0 = triangles + _ids[ind];
                    const uint32_t v00 = (*tria0)[0];
                    const uint32_t v01 = (*tria0)[1];
                    const uint32_t v02 = (*tria0)[2];
                    const vec3<uint32_t>* tria1 = triangles + _ids[++ind];
                    const uint32_t v10 = (*tria1)[0];
                    const uint32_t v11 = (*tria1)[1];
                    const uint32_t v12 = (*tria1)[2];
                    const vec3<uint32_t>* tria2 = triangles + _ids[++ind];
                    const uint32_t v20 = (*tria2)[0];
                    const uint32_t v21 = (*tria2)[1];
                    const uint32_t v22 = (*tria2)[2];
                    const vec3<uint32_t>* tria3 = triangles + _ids[++ind];
                    const uint32_t v30 = (*tria3)[0];
                    const uint32_t v31 = (*tria3)[1];
                    const uint32_t v32 = (*tria3)[2];
                    ++ind;

                    _mm_prefetch((char *)(vertices + v00), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v10), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v20), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v30), _MM_HINT_T0);

                    _mm_prefetch((char *)(vertices + v01), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v11), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v21), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v31), _MM_HINT_T0);

                    _mm_prefetch((char *)(vertices + v02), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v12), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v22), _MM_HINT_T0);
                    _mm_prefetch((char *)(vertices + v32), _MM_HINT_T0);

                    acc.v0[0] = float4(vertices[v00][0], vertices[v10][0], vertices[v20][0], vertices[v30][0]);
                    acc.v0[1] = float4(vertices[v00][1], vertices[v10][1], vertices[v20][1], vertices[v30][1]);
                    acc.v0[2] = float4(vertices[v00][2], vertices[v10][2], vertices[v20][2], vertices[v30][2]);


                    acc.v1[0] = float4(vertices[v01][0], vertices[v11][0], vertices[v21][0], vertices[v31][0]);
                    acc.v1[1] = float4(vertices[v01][1], vertices[v11][1], vertices[v21][1], vertices[v31][1]);
                    acc.v1[2] = float4(vertices[v01][2], vertices[v11][2], vertices[v21][2], vertices[v31][2]);


                    acc.v2[0] = float4(vertices[v02][0], vertices[v12][0], vertices[v22][0], vertices[v32][0]);
                    acc.v2[1] = float4(vertices[v02][1], vertices[v12][1], vertices[v22][1], vertices[v32][1]);
                    acc.v2[2] = float4(vertices[v02][2], vertices[v12][2], vertices[v22][2], vertices[v32][2]);

                    hit4 local_h;
                    intersect_woop(acc, woop_precomputation, ray_origin, t_near, t_far, local_h);

                    if (any(local_h.found))
                      {
                      for (int m = 0; m < 4; ++m)
                        {
                        if (local_h.found[m])
                          {
                          const uint32_t tria_offset = k * 4 + m;
                          if (tria_offset < node.nr_of_primitives[i])
                            {
                            const uint32_t triangle_id = _ids[index[i] + tria_offset];
                            hit h;
                            h.found = -1;
                            h.distance = local_h.distance[m];
                            h.u = local_h.u[m];
                            h.v = local_h.v[m];
                            hits.push_back(h);
                            triangle_ids.push_back(triangle_id);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              else
                stack.push_back((int32_t)node.child[i]);
              }
            }
          }
        else
          {
          for (int i = 0; i < 4; ++i)
            {
            if (isct[i])
              {
              stack.push_back((int32_t)node.child[i]);
              }
            }
          }
        }
      }
    return hits;
    }

  JTKQBVHINLINE hit qbvh::find_closest_triangle(uint32_t& triangle_id, const vec3<float>& point, const vec3<uint32_t>* triangles, const vec3<float>* vertices) const
    {
    hit h;
    h.found = 0;
    h.distance = std::numeric_limits<float>::max();

    vec3<float4> pt(point[0], point[1], point[2]);

    std::vector<int32_t> stack;
    stack.reserve(512);
    stack.push_back((int32_t)root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();

      auto distance_to_node = distance_sqr(node.bbox, pt);
      uint32_t local_order[4] = { 0, 1, 2, 3 };
      std::sort(std::begin(local_order), std::end(local_order), [&](uint32_t a, uint32_t b)
        {
        return distance_to_node[a] > distance_to_node[b];
        });

      const int4 is_leaf = node.child & sign_mask;
      const int contains_at_least_one_leaf = any(is_leaf);
      if (contains_at_least_one_leaf)
        {
        const int4 index = node.child & bit_mask;
        for (int i = 0; i < 4; ++i)
          {
          const uint32_t o = local_order[i];
          if (is_leaf[o])
            {
            const uint32_t d = ((uint32_t)node.nr_of_primitives[o] + 3) >> 2;
            if (d) // not empty
              {
              uint64_t ind = index[o];
              for (uint16_t k = 0; k < d; ++k)
                {
                woop_triangle acc;
                const vec3<uint32_t>* tria0 = triangles + _ids[ind];
                const uint32_t v00 = (*tria0)[0];
                const uint32_t v01 = (*tria0)[1];
                const uint32_t v02 = (*tria0)[2];
                const vec3<uint32_t>* tria1 = triangles + _ids[++ind];
                const uint32_t v10 = (*tria1)[0];
                const uint32_t v11 = (*tria1)[1];
                const uint32_t v12 = (*tria1)[2];
                const vec3<uint32_t>* tria2 = triangles + _ids[++ind];
                const uint32_t v20 = (*tria2)[0];
                const uint32_t v21 = (*tria2)[1];
                const uint32_t v22 = (*tria2)[2];
                const vec3<uint32_t>* tria3 = triangles + _ids[++ind];
                const uint32_t v30 = (*tria3)[0];
                const uint32_t v31 = (*tria3)[1];
                const uint32_t v32 = (*tria3)[2];
                ++ind;

                _mm_prefetch((char *)(vertices + v00), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v10), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v20), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v30), _MM_HINT_T0);

                _mm_prefetch((char *)(vertices + v01), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v11), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v21), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v31), _MM_HINT_T0);

                _mm_prefetch((char *)(vertices + v02), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v12), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v22), _MM_HINT_T0);
                _mm_prefetch((char *)(vertices + v32), _MM_HINT_T0);

                acc.v0[0] = float4(vertices[v00][0], vertices[v10][0], vertices[v20][0], vertices[v30][0]);
                acc.v0[1] = float4(vertices[v00][1], vertices[v10][1], vertices[v20][1], vertices[v30][1]);
                acc.v0[2] = float4(vertices[v00][2], vertices[v10][2], vertices[v20][2], vertices[v30][2]);


                acc.v1[0] = float4(vertices[v01][0], vertices[v11][0], vertices[v21][0], vertices[v31][0]);
                acc.v1[1] = float4(vertices[v01][1], vertices[v11][1], vertices[v21][1], vertices[v31][1]);
                acc.v1[2] = float4(vertices[v01][2], vertices[v11][2], vertices[v21][2], vertices[v31][2]);


                acc.v2[0] = float4(vertices[v02][0], vertices[v12][0], vertices[v22][0], vertices[v32][0]);
                acc.v2[1] = float4(vertices[v02][1], vertices[v12][1], vertices[v22][1], vertices[v32][1]);
                acc.v2[2] = float4(vertices[v02][2], vertices[v12][2], vertices[v22][2], vertices[v32][2]);

                distance4 local_dist;
                distance_sqr(acc, pt, local_dist);

                const bool4 mask = (local_dist.distance_sqr < float4(h.distance));
                if (any(mask))
                  {
                  for (int m = 0; m < 4; ++m)
                    {
                    if (mask[m] && (local_dist.distance_sqr[m] < h.distance))
                      {
                      h.found = -1;
                      h.distance = local_dist.distance_sqr[m];
                      triangle_id = _ids[index[o] + k * 4 + m];
                      h.u = local_dist.u[m];
                      h.v = local_dist.v[m];
                      }
                    }
                  }
                }
              }
            }
          else
            {
            if (distance_to_node[o] < h.distance)
              stack.push_back((int32_t)node.child[o]);
            }
          }
        }

      else
        {
        for (int i = 0; i < 4; ++i)
          {
          const uint32_t o = local_order[i];
          if (distance_to_node[o] < h.distance)
            stack.push_back((int32_t)node.child[o]);
          }
        }
      }
    h.distance = std::sqrt(h.distance);
    return h;
    }


  JTKQBVHINLINE void qbvh::_build(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb)
    {
    assert(props.leaf_size < uint32_t(std::numeric_limits<uint16_t>::max() - 3));
    // max number of nodes is 4n-3 with n the number of triangles
    uint32_t sz;
    //nodes.reserve(4 * nr_of_items / props.leaf_size);

    std::vector<uint32_t> local_ids;
    local_ids.reserve(nr_of_items + 3);
    local_ids.resize(nr_of_items);
    std::iota(local_ids.begin(), local_ids.end(), 0);

    start = local_ids.begin();
    root_id = construct_tree(sz, voxels, total_bb, centroid_bb, local_ids.begin(), local_ids.end());
    if ((int32_t)root_id < 0)
      {
      qbvh_node n;
      n.axis0 = 0;
      n.axis1 = 0;
      n.axis2 = 0;
      n.nr_of_primitives[0] = (uint16_t)sz;
      n.nr_of_primitives[1] = 0;
      n.nr_of_primitives[2] = 0;
      n.nr_of_primitives[3] = 0;
      n.bbox[0] = float4(1);
      n.bbox[1] = float4(1);
      n.bbox[2] = float4(1);
      n.bbox[3] = float4(-1);
      n.bbox[4] = float4(-1);
      n.bbox[5] = float4(-1);

      const auto it = local_ids.begin();
      const auto it_end = it + sz;
      if (voxels)
        get_bbox(n.bbox, voxels, it, it_end, 0);

      n.child[0] = root_id;
      n.child[1] = std::numeric_limits<int32_t>::min();
      n.child[2] = std::numeric_limits<int32_t>::min();
      n.child[3] = std::numeric_limits<int32_t>::min();
      root_id = 0;
      nodes.push_back(n);
      }

    if (!local_ids.empty())
      {
      auto last_id = local_ids.back();
      local_ids.push_back(last_id);
      local_ids.push_back(last_id);
      local_ids.push_back(last_id);
      }

    if (_ids.empty())
      {
      _ids.swap(local_ids);
      }
    else
      {
      for (auto& li : local_ids)
        {
        const uint32_t v = _ids[li];
        li = v;
        }
      _ids.swap(local_ids);
      }
    }

  JTKQBVHINLINE typename int4::value_type qbvh::construct_tree_single(aligned_vector<qbvh_node>& local_nodes, uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    assert(first <= last);
    if (first == last)
      {
      sz = 0;
      return std::numeric_limits<int32_t>::min();
      }
    sz = (uint32_t)std::distance(first, last);
    if (sz <= props.leaf_size)
      {
      typename int4::value_type index = (typename int4::value_type)std::distance(start, first);
      index |= sign_bit;
      return (typename int4::value_type)index;
      }
    else
      {
      uint8_t dim;
      qbvh_node n;
      std::vector<uint32_t>::iterator mid0, mid1, mid2;

      qbvh_voxel bb_left, bb_right, centroid_left, centroid_right;

      qbvh_voxel bb_child[4], centroid_child[4];

      sah_optimized<16>(mid0, bb_left, bb_right, centroid_left, centroid_right, dim, total_bb, centroid_bb, voxels, first, last);
      n.axis0 = dim;
      sah_optimized<16>(mid1, bb_child[0], bb_child[1], centroid_child[0], centroid_child[1], dim, bb_left, centroid_left, voxels, first, mid0);
      n.axis1 = dim;
      sah_optimized<16>(mid2, bb_child[2], bb_child[3], centroid_child[2], centroid_child[3], dim, bb_right, centroid_right, voxels, mid0, last);
      n.axis2 = dim;

      uint32_t sizes[4];

      n.child[0] = construct_tree_single(local_nodes, sizes[0], voxels, bb_child[0], centroid_child[0], first, mid1);
      n.child[1] = construct_tree_single(local_nodes, sizes[1], voxels, bb_child[1], centroid_child[1], mid1, mid0);
      n.child[2] = construct_tree_single(local_nodes, sizes[2], voxels, bb_child[2], centroid_child[2], mid0, mid2);
      n.child[3] = construct_tree_single(local_nodes, sizes[3], voxels, bb_child[3], centroid_child[3], mid2, last);

      for (int i = 0; i < 4; ++i)
        {
        if (n.child[i] >= 0)
          {
          unite_four_aabbs(n.bbox, local_nodes[n.child[i]].bbox, i);
          n.nr_of_primitives[i] = 0;
          }
        else
          {
          n.nr_of_primitives[i] = (uint16_t)sizes[i];
          n.bbox[0][i] = bb_child[i].bbox_min[0];
          n.bbox[1][i] = bb_child[i].bbox_min[1];
          n.bbox[2][i] = bb_child[i].bbox_min[2];
          n.bbox[3][i] = bb_child[i].bbox_max[0];
          n.bbox[4][i] = bb_child[i].bbox_max[1];
          n.bbox[5][i] = bb_child[i].bbox_max[2];
          }
        }
      uint32_t node_id = (uint32_t)local_nodes.size();
      local_nodes.push_back(n);
      return (int32_t)node_id;
      }
    }


  JTKQBVHINLINE typename int4::value_type qbvh::construct_tree_prep(std::vector<tree_stack>& stack, uint32_t& sz, uint32_t& stack_index, uint32_t& node_index, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last, int level, bool prep)
    {
    assert(first <= last);
    if (first == last)
      {
      sz = 0;
      return std::numeric_limits<int32_t>::min();
      }
    sz = (uint32_t)std::distance(first, last);
    if (sz <= props.leaf_size)
      {
      typename int4::value_type index = (typename int4::value_type)std::distance(start, first);
      index |= sign_bit;
      return (typename int4::value_type)index;
      }
    else
      {
      tree_stack st;
      uint32_t stack_id;


      qbvh_node n;
      uint32_t sizes[4];

      if (prep)
        {
        uint8_t dim;
        stack_id = (uint32_t)stack.size();
        stack.emplace_back();

        qbvh_voxel bb_left, bb_right, centroid_left, centroid_right;
        st.first = first;
        st.last = last;
        sah_parallel<16>(st.mid0, bb_left, bb_right, centroid_left, centroid_right, dim, total_bb, centroid_bb, voxels, st.first, st.last);
        n.axis0 = dim;
        sah_parallel<16>(st.mid1, st.bb_child[0], st.bb_child[1], st.centroid_child[0], st.centroid_child[1], dim, bb_left, centroid_left, voxels, st.first, st.mid0);
        n.axis1 = dim;
        sah_parallel<16>(st.mid2, st.bb_child[2], st.bb_child[3], st.centroid_child[2], st.centroid_child[3], dim, bb_right, centroid_right, voxels, st.mid0, st.last);
        n.axis2 = dim;
        st.level = level;
        }
      else
        {
        stack_id = stack_index++;
        st = stack[stack_id];
        }
      assert(st.level == level);
      if (level == props.parallel_level)
        {
        if (!prep)
          {
          for (int i = 0; i < 4; ++i)
            {
            n.child[i] = st.child[i];
            sizes[i] = st.sizes[i];
            }
          }
        }
      else
        {
        n.child[0] = construct_tree_prep(stack, sizes[0], stack_index, node_index, voxels, st.bb_child[0], st.centroid_child[0], st.first, st.mid1, level + 1, prep);
        n.child[1] = construct_tree_prep(stack, sizes[1], stack_index, node_index, voxels, st.bb_child[1], st.centroid_child[1], st.mid1, st.mid0, level + 1, prep);
        n.child[2] = construct_tree_prep(stack, sizes[2], stack_index, node_index, voxels, st.bb_child[2], st.centroid_child[2], st.mid0, st.mid2, level + 1, prep);
        n.child[3] = construct_tree_prep(stack, sizes[3], stack_index, node_index, voxels, st.bb_child[3], st.centroid_child[3], st.mid2, st.last, level + 1, prep);
        st.child[0] = n.child[0];
        st.child[1] = n.child[1];
        st.child[2] = n.child[2];
        st.child[3] = n.child[3];
        st.sizes[0] = sizes[0];
        st.sizes[1] = sizes[1];
        st.sizes[2] = sizes[2];
        st.sizes[3] = sizes[3];
        }
      if (!prep)
        {
        for (int i = 0; i < 4; ++i)
          {
          if (n.child[i] >= 0)
            {
            unite_four_aabbs(n.bbox, nodes[n.child[i]].bbox, i);
            }
          else
            {
            n.nr_of_primitives[i] = uint16_t(sizes[i]);
            n.bbox[0][i] = st.bb_child[i].bbox_min[0];
            n.bbox[1][i] = st.bb_child[i].bbox_min[1];
            n.bbox[2][i] = st.bb_child[i].bbox_min[2];
            n.bbox[3][i] = st.bb_child[i].bbox_max[0];
            n.bbox[4][i] = st.bb_child[i].bbox_max[1];
            n.bbox[5][i] = st.bb_child[i].bbox_max[2];
            }
          }
        }
      if (prep)
        {
        uint32_t node_id = (uint32_t)nodes.size();
        nodes.push_back(n);
        stack[stack_id] = st;
        return node_id;
        }
      else
        {
        n.axis0 = nodes[node_index].axis0;
        n.axis1 = nodes[node_index].axis1;
        n.axis2 = nodes[node_index].axis2;
        nodes[node_index] = n;
        ++node_index;
        return node_index - 1;
        }
      //return 0; // unreachable
      }
    }


  JTKQBVHINLINE typename int4::value_type qbvh::construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    std::vector<tree_stack> stack;
    uint32_t stack_index = 0;
    uint32_t node_index = 0;
    const auto id = construct_tree_prep(stack, sz, stack_index, node_index, voxels, total_bb, centroid_bb, first, last, 0, true);
    const uint32_t stack_size = (uint32_t)stack.size();

    std::vector<uint32_t> subtrees;
    subtrees.reserve(stack_size);
    for (uint32_t i = 0; i < stack_size; ++i)
      {
      if (stack[i].level == props.parallel_level)
        subtrees.push_back(i);
      }

    std::sort(subtrees.begin(), subtrees.end(), [&](uint32_t left, uint32_t right)
      {
      return std::distance(stack[left].first, stack[left].last) > std::distance(stack[right].first, stack[right].last);
      });

    const uint32_t subtree_size = (uint32_t)subtrees.size();

    std::vector<aligned_vector<qbvh_node>> local_nodes(subtree_size * 4, nodes);
    int32_t front_size = (int32_t)nodes.size();
    parallel_for(uint32_t(0), subtree_size * 4, [&](uint32_t i)
      {
      const uint32_t subtree_id = i / 4;
      const uint32_t j = i & 3;
      tree_stack& st = stack[subtrees[subtree_id]];
      switch (j)
        {
        case 0: st.child[0] = construct_tree_single(local_nodes[i], st.sizes[0], voxels, st.bb_child[0], st.centroid_child[0], st.first, st.mid1); break;
        case 1: st.child[1] = construct_tree_single(local_nodes[i], st.sizes[1], voxels, st.bb_child[1], st.centroid_child[1], st.mid1, st.mid0); break;
        case 2: st.child[2] = construct_tree_single(local_nodes[i], st.sizes[2], voxels, st.bb_child[2], st.centroid_child[2], st.mid0, st.mid2); break;
        case 3: st.child[3] = construct_tree_single(local_nodes[i], st.sizes[3], voxels, st.bb_child[3], st.centroid_child[3], st.mid2, st.last); break;
        }
      });

    uint64_t reserve_nodes_size = front_size;
    for (uint32_t i = 0; i < subtree_size * 4; ++i)
      {
      reserve_nodes_size += (local_nodes[i].size() - front_size);
      }
    nodes.reserve(reserve_nodes_size);

    for (uint32_t i = 0; i < subtree_size * 4; ++i)
      {
      const uint32_t subtree_id = i / 4;
      const uint32_t j = i & 3;
      tree_stack& st = stack[subtrees[subtree_id]];
      int32_t offset = (int32_t)nodes.size() - front_size;

      switch (j)
        {
        case 0: if (st.child[0] >= 0) st.child[0] += offset; break;
        case 1: if (st.child[1] >= 0) st.child[1] += offset; break;
        case 2: if (st.child[2] >= 0) st.child[2] += offset; break;
        case 3: if (st.child[3] >= 0) st.child[3] += offset; break;
        }

      for (int32_t k = front_size; k < local_nodes[i].size(); ++k)
        {
        auto& n = local_nodes[i][k];
        if (n.child[0] >= front_size)
          n.child[0] += offset;
        if (n.child[1] >= front_size)
          n.child[1] += offset;
        if (n.child[2] >= front_size)
          n.child[2] += offset;
        if (n.child[3] >= front_size)
          n.child[3] += offset;
        nodes.push_back(n);
        }
      }

    const auto id2 = construct_tree_prep(stack, sz, stack_index, node_index, voxels, total_bb, centroid_bb, first, last, 0, false);
    (void)id2; // avoid unused variable warning in release build
    assert(id == id2);

    return id;
    }

  /////////////////////////////////////////////////////////////////////////
  // sphere_qbvh
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHINLINE sphere_qbvh::sphere_qbvh(const vec3<float>* origins, const float* radii, uint32_t nr_of_spheres)
    {
    qbvh_voxel total_bb, centroid_bb;
    auto qbvh_voxels = build_sphere_qbvh_voxels(total_bb, centroid_bb, origins, radii, nr_of_spheres);
    _build(qbvh_voxels, nr_of_spheres, total_bb, centroid_bb);
    delete[] qbvh_voxels;
    }

  JTKQBVHINLINE sphere_qbvh::sphere_qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const std::vector<uint32_t>& ids, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr) : props(pr)
    {
    assert(props.leaf_size > 1);
    _ids = ids;
    _build(voxels, nr_of_items, total_bb, centroid_bb);
    }

  JTKQBVHINLINE sphere_qbvh::sphere_qbvh(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb, properties pr) : props(pr)
    {
    assert(props.leaf_size > 1);
    _build(voxels, nr_of_items, total_bb, centroid_bb);
    }

  JTKQBVHINLINE spherehit sphere_qbvh::find_closest_sphere(uint32_t& sphere_id, ray r, const vec3<float>* origins, const float* radii) const
    {
    spherehit h;
    h.found = 0;
    h.distance = std::numeric_limits<float>::max();

    float r_dir_length = std::sqrt(r.dir[0] * r.dir[0] + r.dir[1] * r.dir[1] + r.dir[2] * r.dir[2]);
    r.dir = r.dir / r_dir_length;

    vec3<float4> ray_origin(r.orig[0], r.orig[1], r.orig[2]);
    vec3<float4> ray_dir(r.dir[0], r.dir[1], r.dir[2]);

    int32_t ray_dir_sign[3];
    float4 ray_inverse_dir_tmp = reciprocal(r.dir);
    if (fabs(ray_inverse_dir_tmp[0]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[0] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[1]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[1] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[2]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[2] = std::numeric_limits<float>::max();
    vec3<float4> ray_inverse_dir(ray_inverse_dir_tmp[0], ray_inverse_dir_tmp[1], ray_inverse_dir_tmp[2]);
    float4 t_near(r.t_near);
    float4 t_far(r.t_far);
    ray_dir_sign[0] = r.dir[0] < 0 ? 1 : 0;
    ray_dir_sign[1] = r.dir[1] < 0 ? 1 : 0;
    ray_dir_sign[2] = r.dir[2] < 0 ? 1 : 0;

    std::vector<int32_t> stack;
    stack.reserve(512);
    stack.push_back((int32_t)root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();

      const uint32_t order_pos = (ray_dir_sign[node.axis0] << 2) | (ray_dir_sign[node.axis1] << 1) | (ray_dir_sign[node.axis2]);
      const uint32_t* local_order = order[order_pos];

      auto isct = intersect(node.bbox, ray_origin, t_near, t_far, ray_dir_sign, ray_inverse_dir);

      if (any(isct))
        {
        const int4 is_leaf = node.child & sign_mask;
        const int contains_at_least_one_leaf = any(is_leaf);
        if (contains_at_least_one_leaf)
          {
          const int4 index = node.child & bit_mask;
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              {
              if (is_leaf[o])
                {
                const uint32_t d = ((uint32_t)node.nr_of_primitives[o] + 3) >> 2;
                if (d) // not empty
                  {
                  uint64_t ind = index[o];
                  for (uint16_t k = 0; k < d; ++k)
                    {

                    const vec3<float>* orig0 = origins + _ids[ind];
                    const float* r0 = radii + _ids[ind];
                    const vec3<float>* orig1 = origins + _ids[++ind];
                    const float* r1 = radii + _ids[ind];
                    const vec3<float>* orig2 = origins + _ids[++ind];
                    const float* r2 = radii + _ids[ind];
                    const vec3<float>* orig3 = origins + _ids[++ind];
                    const float* r3 = radii + _ids[ind];
                    ++ind;

                    const float4 o0((*orig0)[0], (*orig1)[0], (*orig2)[0], (*orig3)[0]);
                    const float4 o1((*orig0)[1], (*orig1)[1], (*orig2)[1], (*orig3)[1]);
                    const float4 o2((*orig0)[2], (*orig1)[2], (*orig2)[2], (*orig3)[2]);

                    vec3<float4> origin(o0, o1, o2);
                    float4 rad(*r0, *r1, *r2, *r3);

                    spherehit4 local_h;
                    intersect_sphere(origin, rad, ray_origin, ray_dir, t_near, t_far, local_h);
                    const bool4 mask = (abs(local_h.distance) < float4(std::abs(h.distance))) & local_h.found;
                    if (any(mask))
                      {
                      for (int m = 0; m < 4; ++m)
                        {
                        if (mask[m])
                          {
                          if (std::abs(local_h.distance[m]) < std::abs(h.distance))
                            {
                            h.found = -1;
                            h.distance = local_h.distance[m];
                            sphere_id = _ids[index[o] + k * 4 + m];
                            (h.distance > 0) ? t_far = float4(h.distance) : t_near = float4(h.distance);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              else
                stack.push_back((int32_t)node.child[o]);
              }
            }
          }
        else
          {
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              {
              stack.push_back((int32_t)node.child[o]);
              }
            }
          }
        }
      }
    return h;
    }

  JTKQBVHINLINE std::vector<spherehit> sphere_qbvh::find_all_spheres(std::vector<uint32_t>& sphere_ids, ray r, const vec3<float>* origins, const float* radii) const
    {
    std::vector<spherehit> hits;

    float r_dir_length = std::sqrt(r.dir[0] * r.dir[0] + r.dir[1] * r.dir[1] + r.dir[2] * r.dir[2]);
    r.dir = r.dir / r_dir_length;

    vec3<float4> ray_origin(r.orig[0], r.orig[1], r.orig[2]);
    vec3<float4> ray_dir(r.dir[0], r.dir[1], r.dir[2]);

    int32_t ray_dir_sign[3];
    float4 ray_inverse_dir_tmp = reciprocal(r.dir);
    if (fabs(ray_inverse_dir_tmp[0]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[0] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[1]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[1] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[2]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[2] = std::numeric_limits<float>::max();
    vec3<float4> ray_inverse_dir(ray_inverse_dir_tmp[0], ray_inverse_dir_tmp[1], ray_inverse_dir_tmp[2]);
    float4 t_near(r.t_near);
    float4 t_far(r.t_far);
    ray_dir_sign[0] = r.dir[0] < 0 ? 1 : 0;
    ray_dir_sign[1] = r.dir[1] < 0 ? 1 : 0;
    ray_dir_sign[2] = r.dir[2] < 0 ? 1 : 0;

    std::vector<int32_t> stack;
    stack.reserve(512);
    stack.push_back((int32_t)root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();

      auto isct = intersect(node.bbox, ray_origin, t_near, t_far, ray_dir_sign, ray_inverse_dir);

      if (any(isct))
        {
        const int4 is_leaf = node.child & sign_mask;
        const int contains_at_least_one_leaf = any(is_leaf);
        if (contains_at_least_one_leaf)
          {
          const int4 index = node.child & bit_mask;
          for (int i = 0; i < 4; ++i)
            {
            if (isct[i])
              {
              if (is_leaf[i])
                {
                const uint32_t d = ((uint32_t)node.nr_of_primitives[i] + 3) >> 2;
                if (d) // not empty
                  {
                  uint64_t ind = index[i];
                  for (uint16_t k = 0; k < d; ++k)
                    {

                    const vec3<float>* orig0 = origins + _ids[ind];
                    const float* r0 = radii + _ids[ind];
                    const vec3<float>* orig1 = origins + _ids[++ind];
                    const float* r1 = radii + _ids[ind];
                    const vec3<float>* orig2 = origins + _ids[++ind];
                    const float* r2 = radii + _ids[ind];
                    const vec3<float>* orig3 = origins + _ids[++ind];
                    const float* r3 = radii + _ids[ind];
                    ++ind;

                    const float4 o0((*orig0)[0], (*orig1)[0], (*orig2)[0], (*orig3)[0]);
                    const float4 o1((*orig0)[1], (*orig1)[1], (*orig2)[1], (*orig3)[1]);
                    const float4 o2((*orig0)[2], (*orig1)[2], (*orig2)[2], (*orig3)[2]);

                    vec3<float4> origin(o0, o1, o2);
                    float4 rad(*r0, *r1, *r2, *r3);

                    spherehit4 local_h;
                    intersect_sphere(origin, rad, ray_origin, ray_dir, t_near, t_far, local_h);
                    if (any(local_h.found))
                      {
                      for (int m = 0; m < 4; ++m)
                        {
                        if (local_h.found[m])
                          {
                          const uint32_t off = k * 4 + m;
                          if (off < node.nr_of_primitives[i])
                            {
                            const uint32_t sphere_id = _ids[index[i] + off];
                            spherehit h;
                            h.found = -1;
                            h.distance = local_h.distance[m];
                            hits.push_back(h);
                            sphere_ids.push_back(sphere_id);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              else
                stack.push_back((int32_t)node.child[i]);
              }
            }
          }
        else
          {
          for (int i = 0; i < 4; ++i)
            {
            if (isct[i])
              {
              stack.push_back((int32_t)node.child[i]);
              }
            }
          }
        }
      }
    return hits;
    }

  JTKQBVHINLINE void sphere_qbvh::_build(const qbvh_voxel* voxels, uint32_t nr_of_items, const qbvh_voxel& total_bb, const qbvh_voxel& centroid_bb)
    {
    assert(props.leaf_size < uint32_t(std::numeric_limits<uint16_t>::max() - 3));
    // max number of nodes is 4n-3 with n the number of triangles
    uint32_t sz;
    //nodes.reserve(4 * nr_of_items / props.leaf_size);

    std::vector<uint32_t> local_ids;
    local_ids.reserve(nr_of_items + 3);
    local_ids.resize(nr_of_items);
    std::iota(local_ids.begin(), local_ids.end(), 0);

    start = local_ids.begin();
    root_id = construct_tree(sz, voxels, total_bb, centroid_bb, local_ids.begin(), local_ids.end());
    if ((int32_t)root_id < 0)
      {
      qbvh_node n;
      n.axis0 = 0;
      n.axis1 = 0;
      n.axis2 = 0;
      n.nr_of_primitives[0] = (uint16_t)sz;
      n.nr_of_primitives[1] = 0;
      n.nr_of_primitives[2] = 0;
      n.nr_of_primitives[3] = 0;
      n.bbox[0] = float4(1);
      n.bbox[1] = float4(1);
      n.bbox[2] = float4(1);
      n.bbox[3] = float4(-1);
      n.bbox[4] = float4(-1);
      n.bbox[5] = float4(-1);

      const auto it = local_ids.begin();
      const auto it_end = it + sz;
      if (voxels)
        get_bbox(n.bbox, voxels, it, it_end, 0);

      n.child[0] = root_id;
      n.child[1] = std::numeric_limits<int32_t>::min();
      n.child[2] = std::numeric_limits<int32_t>::min();
      n.child[3] = std::numeric_limits<int32_t>::min();
      root_id = 0;
      nodes.push_back(n);
      }

    if (!local_ids.empty())
      {
      auto last_id = local_ids.back();
      local_ids.push_back(last_id);
      local_ids.push_back(last_id);
      local_ids.push_back(last_id);
      }

    if (_ids.empty())
      {
      _ids.swap(local_ids);
      }
    else
      {
      for (auto& li : local_ids)
        {
        const uint32_t v = _ids[li];
        li = v;
        }
      _ids.swap(local_ids);
      }
    }

  JTKQBVHINLINE typename int4::value_type sphere_qbvh::construct_tree_single(aligned_vector<qbvh_node>& local_nodes, uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    assert(first <= last);
    if (first == last)
      {
      sz = 0;
      return std::numeric_limits<int32_t>::min();
      }
    sz = (uint32_t)std::distance(first, last);
    if (sz <= props.leaf_size)
      {
      typename int4::value_type index = (typename int4::value_type)std::distance(start, first);
      index |= sign_bit;
      return (typename int4::value_type)index;
      }
    else
      {
      uint8_t dim;
      qbvh_node n;
      std::vector<uint32_t>::iterator mid0, mid1, mid2;

      qbvh_voxel bb_left, bb_right, centroid_left, centroid_right;

      qbvh_voxel bb_child[4], centroid_child[4];

      sah_optimized<16>(mid0, bb_left, bb_right, centroid_left, centroid_right, dim, total_bb, centroid_bb, voxels, first, last);
      n.axis0 = dim;
      sah_optimized<16>(mid1, bb_child[0], bb_child[1], centroid_child[0], centroid_child[1], dim, bb_left, centroid_left, voxels, first, mid0);
      n.axis1 = dim;
      sah_optimized<16>(mid2, bb_child[2], bb_child[3], centroid_child[2], centroid_child[3], dim, bb_right, centroid_right, voxels, mid0, last);
      n.axis2 = dim;

      uint32_t sizes[4];

      n.child[0] = construct_tree_single(local_nodes, sizes[0], voxels, bb_child[0], centroid_child[0], first, mid1);
      n.child[1] = construct_tree_single(local_nodes, sizes[1], voxels, bb_child[1], centroid_child[1], mid1, mid0);
      n.child[2] = construct_tree_single(local_nodes, sizes[2], voxels, bb_child[2], centroid_child[2], mid0, mid2);
      n.child[3] = construct_tree_single(local_nodes, sizes[3], voxels, bb_child[3], centroid_child[3], mid2, last);

      for (int i = 0; i < 4; ++i)
        {
        if (n.child[i] >= 0)
          {
          unite_four_aabbs(n.bbox, local_nodes[n.child[i]].bbox, i);
          n.nr_of_primitives[i] = 0;
          }
        else
          {
          n.nr_of_primitives[i] = (uint16_t)sizes[i];
          n.bbox[0][i] = bb_child[i].bbox_min[0];
          n.bbox[1][i] = bb_child[i].bbox_min[1];
          n.bbox[2][i] = bb_child[i].bbox_min[2];
          n.bbox[3][i] = bb_child[i].bbox_max[0];
          n.bbox[4][i] = bb_child[i].bbox_max[1];
          n.bbox[5][i] = bb_child[i].bbox_max[2];
          }
        }
      uint32_t node_id = (uint32_t)local_nodes.size();
      local_nodes.push_back(n);
      return (int32_t)node_id;
      }
    }


  JTKQBVHINLINE typename int4::value_type sphere_qbvh::construct_tree_prep(std::vector<tree_stack>& stack, uint32_t& sz, uint32_t& stack_index, uint32_t& node_index, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last, int level, bool prep)
    {
    assert(first <= last);
    if (first == last)
      {
      sz = 0;
      return std::numeric_limits<int32_t>::min();
      }
    sz = (uint32_t)std::distance(first, last);
    if (sz <= props.leaf_size)
      {
      typename int4::value_type index = (typename int4::value_type)std::distance(start, first);
      index |= sign_bit;
      return (typename int4::value_type)index;
      }
    else
      {
      tree_stack st;
      uint32_t stack_id;


      qbvh_node n;
      uint32_t sizes[4];

      if (prep)
        {
        uint8_t dim;
        stack_id = (uint32_t)stack.size();
        stack.emplace_back();

        qbvh_voxel bb_left, bb_right, centroid_left, centroid_right;
        st.first = first;
        st.last = last;
        sah_parallel<16>(st.mid0, bb_left, bb_right, centroid_left, centroid_right, dim, total_bb, centroid_bb, voxels, st.first, st.last);
        n.axis0 = dim;
        sah_parallel<16>(st.mid1, st.bb_child[0], st.bb_child[1], st.centroid_child[0], st.centroid_child[1], dim, bb_left, centroid_left, voxels, st.first, st.mid0);
        n.axis1 = dim;
        sah_parallel<16>(st.mid2, st.bb_child[2], st.bb_child[3], st.centroid_child[2], st.centroid_child[3], dim, bb_right, centroid_right, voxels, st.mid0, st.last);
        n.axis2 = dim;
        st.level = level;
        }
      else
        {
        stack_id = stack_index++;
        st = stack[stack_id];
        }
      assert(st.level == level);
      if (level == props.parallel_level)
        {
        if (!prep)
          {
          for (int i = 0; i < 4; ++i)
            {
            n.child[i] = st.child[i];
            sizes[i] = st.sizes[i];
            }
          }
        }
      else
        {
        n.child[0] = construct_tree_prep(stack, sizes[0], stack_index, node_index, voxels, st.bb_child[0], st.centroid_child[0], st.first, st.mid1, level + 1, prep);
        n.child[1] = construct_tree_prep(stack, sizes[1], stack_index, node_index, voxels, st.bb_child[1], st.centroid_child[1], st.mid1, st.mid0, level + 1, prep);
        n.child[2] = construct_tree_prep(stack, sizes[2], stack_index, node_index, voxels, st.bb_child[2], st.centroid_child[2], st.mid0, st.mid2, level + 1, prep);
        n.child[3] = construct_tree_prep(stack, sizes[3], stack_index, node_index, voxels, st.bb_child[3], st.centroid_child[3], st.mid2, st.last, level + 1, prep);
        st.child[0] = n.child[0];
        st.child[1] = n.child[1];
        st.child[2] = n.child[2];
        st.child[3] = n.child[3];
        st.sizes[0] = sizes[0];
        st.sizes[1] = sizes[1];
        st.sizes[2] = sizes[2];
        st.sizes[3] = sizes[3];
        }
      if (!prep)
        {
        for (int i = 0; i < 4; ++i)
          {
          if (n.child[i] >= 0)
            {
            unite_four_aabbs(n.bbox, nodes[n.child[i]].bbox, i);
            }
          else
            {
            n.nr_of_primitives[i] = uint16_t(sizes[i]);
            n.bbox[0][i] = st.bb_child[i].bbox_min[0];
            n.bbox[1][i] = st.bb_child[i].bbox_min[1];
            n.bbox[2][i] = st.bb_child[i].bbox_min[2];
            n.bbox[3][i] = st.bb_child[i].bbox_max[0];
            n.bbox[4][i] = st.bb_child[i].bbox_max[1];
            n.bbox[5][i] = st.bb_child[i].bbox_max[2];
            }
          }
        }
      if (prep)
        {
        uint32_t node_id = (uint32_t)nodes.size();
        nodes.push_back(n);
        stack[stack_id] = st;
        return node_id;
        }
      else
        {
        n.axis0 = nodes[node_index].axis0;
        n.axis1 = nodes[node_index].axis1;
        n.axis2 = nodes[node_index].axis2;
        nodes[node_index] = n;
        ++node_index;
        return node_index - 1;
        }
      //return 0; // unreachable
      }
    }


  JTKQBVHINLINE typename int4::value_type sphere_qbvh::construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    std::vector<tree_stack> stack;
    uint32_t stack_index = 0;
    uint32_t node_index = 0;
    const auto id = construct_tree_prep(stack, sz, stack_index, node_index, voxels, total_bb, centroid_bb, first, last, 0, true);
    const uint32_t stack_size = (uint32_t)stack.size();

    std::vector<uint32_t> subtrees;
    subtrees.reserve(stack_size);
    for (uint32_t i = 0; i < stack_size; ++i)
      {
      if (stack[i].level == props.parallel_level)
        subtrees.push_back(i);
      }

    std::sort(subtrees.begin(), subtrees.end(), [&](uint32_t left, uint32_t right)
      {
      return std::distance(stack[left].first, stack[left].last) > std::distance(stack[right].first, stack[right].last);
      });

    const uint32_t subtree_size = (uint32_t)subtrees.size();

    std::vector<aligned_vector<qbvh_node>> local_nodes(subtree_size * 4, nodes);
    int32_t front_size = (int32_t)nodes.size();
    parallel_for(uint32_t(0), subtree_size * 4, [&](uint32_t i)
      {
      const uint32_t subtree_id = i / 4;
      const uint32_t j = i & 3;
      tree_stack& st = stack[subtrees[subtree_id]];
      switch (j)
        {
        case 0: st.child[0] = construct_tree_single(local_nodes[i], st.sizes[0], voxels, st.bb_child[0], st.centroid_child[0], st.first, st.mid1); break;
        case 1: st.child[1] = construct_tree_single(local_nodes[i], st.sizes[1], voxels, st.bb_child[1], st.centroid_child[1], st.mid1, st.mid0); break;
        case 2: st.child[2] = construct_tree_single(local_nodes[i], st.sizes[2], voxels, st.bb_child[2], st.centroid_child[2], st.mid0, st.mid2); break;
        case 3: st.child[3] = construct_tree_single(local_nodes[i], st.sizes[3], voxels, st.bb_child[3], st.centroid_child[3], st.mid2, st.last); break;
        }
      });

    uint64_t reserve_nodes_size = front_size;
    for (uint32_t i = 0; i < subtree_size * 4; ++i)
      {
      reserve_nodes_size += (local_nodes[i].size() - front_size);
      }
    nodes.reserve(reserve_nodes_size);

    for (uint32_t i = 0; i < subtree_size * 4; ++i)
      {
      const uint32_t subtree_id = i / 4;
      const uint32_t j = i & 3;
      tree_stack& st = stack[subtrees[subtree_id]];
      int32_t offset = (int32_t)nodes.size() - front_size;

      switch (j)
        {
        case 0: if (st.child[0] >= 0) st.child[0] += offset; break;
        case 1: if (st.child[1] >= 0) st.child[1] += offset; break;
        case 2: if (st.child[2] >= 0) st.child[2] += offset; break;
        case 3: if (st.child[3] >= 0) st.child[3] += offset; break;
        }

      for (int32_t k = front_size; k < local_nodes[i].size(); ++k)
        {
        auto& n = local_nodes[i][k];
        if (n.child[0] >= front_size)
          n.child[0] += offset;
        if (n.child[1] >= front_size)
          n.child[1] += offset;
        if (n.child[2] >= front_size)
          n.child[2] += offset;
        if (n.child[3] >= front_size)
          n.child[3] += offset;
        nodes.push_back(n);
        }
      }

    const auto id2 = construct_tree_prep(stack, sz, stack_index, node_index, voxels, total_bb, centroid_bb, first, last, 0, false);
    (void)id2; // avoid unused variable warning in release build
    assert(id == id2);

    return id;
    }
    
    

  JTKQBVHINLINE qbvh_two_level::qbvh_two_level(const qbvh** objects, uint32_t nr_of_objects)
    {
    nodes.reserve(nr_of_objects);
    qbvh_voxel total_bb, centroid_bb;
    qbvh_voxel* voxels = new qbvh_voxel[nr_of_objects];
    total_bb.bbox_min = std::numeric_limits<float>::max();
    total_bb.bbox_max = -std::numeric_limits<float>::max();
    centroid_bb.bbox_min = total_bb.bbox_min;
    centroid_bb.bbox_max = total_bb.bbox_max;
    for (uint32_t t = 0; t < nr_of_objects; ++t)
      {
      const qbvh* current_bvh = objects[t];
      unite_four_aabbs(voxels[t], current_bvh->nodes[current_bvh->root_id].bbox);
      voxels[t].centroid = (voxels[t].bbox_min + voxels[t].bbox_max)*0.5f;
      total_bb.bbox_min = min(total_bb.bbox_min, voxels[t].bbox_min);
      total_bb.bbox_max = max(total_bb.bbox_max, voxels[t].bbox_max);
      centroid_bb.bbox_min = min(centroid_bb.bbox_min, voxels[t].centroid);
      centroid_bb.bbox_max = max(centroid_bb.bbox_max, voxels[t].centroid);
      }
    object_ids.resize(nr_of_objects);
    std::iota(object_ids.begin(), object_ids.end(), 0);
    uint32_t sz;
    start = object_ids.begin();
    root_id = construct_tree(sz, voxels, total_bb, centroid_bb, object_ids.begin(), object_ids.end());
    if ((int32_t)root_id < 0)
      {
      qbvh_two_level_node n;
      n.axis0 = 0;
      n.axis1 = 0;
      n.axis2 = 0;
      n.bbox[0] = float4(1.f);
      n.bbox[1] = float4(1.f);
      n.bbox[2] = float4(1.f);
      n.bbox[3] = float4(-1.f);
      n.bbox[4] = float4(-1.f);
      n.bbox[5] = float4(-1.f);

      const auto it = object_ids.begin();
      const auto it_end = it + sz;
      get_bbox(n.bbox, voxels, it, it_end, 0);

      n.child[0] = root_id;
      n.child[1] = std::numeric_limits<int32_t>::min();
      n.child[2] = std::numeric_limits<int32_t>::min();
      n.child[3] = std::numeric_limits<int32_t>::min();
      root_id = 0;
      nodes.push_back(n);
      }
    delete[] voxels;
    }

  JTKQBVHINLINE hit qbvh_two_level::find_closest_triangle(uint32_t& triangle_id, uint32_t& object_id, ray r, const qbvh** objects, const vec3<uint32_t>** triangles, const vec3<float>** vertices) const
    {
    hit h;
    h.found = 0;
    h.distance = std::numeric_limits<float>::max();

    vec3<float4> ray_origin(r.orig[0], r.orig[1], r.orig[2]);
    vec3<float4> ray_dir(r.dir[0], r.dir[1], r.dir[2]);

    int32_t ray_dir_sign[3];
    float4 ray_inverse_dir_tmp = reciprocal(r.dir);
    if (fabs(ray_inverse_dir_tmp[0]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[0] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[1]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[1] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[2]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[2] = std::numeric_limits<float>::max();
    vec3<float4> ray_inverse_dir(ray_inverse_dir_tmp[0], ray_inverse_dir_tmp[1], ray_inverse_dir_tmp[2]);
    float4 t_near(r.t_near);
    float4 t_far(r.t_far);
    ray_dir_sign[0] = r.dir[0] < 0 ? 1 : 0;
    ray_dir_sign[1] = r.dir[1] < 0 ? 1 : 0;
    ray_dir_sign[2] = r.dir[2] < 0 ? 1 : 0;

    std::vector<int32_t> stack;
    stack.reserve(32);
    stack.push_back(root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();
      const uint32_t order_pos = (ray_dir_sign[node.axis0] << 2) | (ray_dir_sign[node.axis1] << 1) | (ray_dir_sign[node.axis2]);
      const uint32_t* local_order = order[order_pos];

      auto isct = intersect(node.bbox, ray_origin, t_near, t_far, ray_dir_sign, ray_inverse_dir);
      if (any(isct))
        {
        const int4 is_leaf = node.child & sign_mask;
        const int contains_at_least_one_leaf = any(is_leaf);
        if (contains_at_least_one_leaf)
          {
          const int4 index = node.child & bit_mask;
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              {
              if (is_leaf[o])
                {
                auto object_it = object_ids.begin() + index[o];
                const uint32_t local_object_id = *object_it;
                uint32_t triangle_candidate;
                ray local_r(r);
                local_r.t_far = t_far[o];
                local_r.t_near = t_near[o];
                auto local_hit = objects[local_object_id]->find_closest_triangle(triangle_candidate, local_r, triangles[local_object_id], vertices[local_object_id]);
                if (local_hit.found && std::abs(local_hit.distance) < std::abs(h.distance))
                  {
                  h = local_hit;
                  triangle_id = triangle_candidate;
                  object_id = local_object_id;
                  (h.distance > 0) ? t_far = float4(h.distance) : t_near = float4(h.distance);
                  }
                }
              else
                stack.push_back(node.child[o]);
              }
            }
          }
        else
          {
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              stack.push_back(node.child[o]);
            }
          }
        }
      }

    return h;
    }


  JTKQBVHINLINE int32_t qbvh_two_level::construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    assert(first <= last);
    if (first == last)
      {
      sz = 0;
      return std::numeric_limits<int32_t>::min();
      }
    sz = (uint32_t)std::distance(first, last);
    if (sz <= leaf_size)
      {
      uint32_t index = (uint32_t)std::distance(start, first);
      index |= sign_bit;
      return (int32_t)index;
      }
    else
      {
      uint8_t dim;
      qbvh_two_level_node n;
      std::vector<uint32_t>::iterator mid0, mid1, mid2;

      qbvh_voxel bb_left, bb_right, centroid_left, centroid_right;

      qbvh_voxel bb_child[4], centroid_child[4];

      sah_optimized<4>(mid0, bb_left, bb_right, centroid_left, centroid_right, dim, total_bb, centroid_bb, voxels, first, last);
      n.axis0 = dim;
      sah_optimized<4>(mid1, bb_child[0], bb_child[1], centroid_child[0], centroid_child[1], dim, bb_left, centroid_left, voxels, first, mid0);
      n.axis1 = dim;
      sah_optimized<4>(mid2, bb_child[2], bb_child[3], centroid_child[2], centroid_child[3], dim, bb_right, centroid_right, voxels, mid0, last);
      n.axis2 = dim;

      uint32_t sizes[4];

      n.child[0] = construct_tree(sizes[0], voxels, bb_child[0], centroid_child[0], first, mid1);
      n.child[1] = construct_tree(sizes[1], voxels, bb_child[1], centroid_child[1], mid1, mid0);
      n.child[2] = construct_tree(sizes[2], voxels, bb_child[2], centroid_child[2], mid0, mid2);
      n.child[3] = construct_tree(sizes[3], voxels, bb_child[3], centroid_child[3], mid2, last);

      for (int i = 0; i < 4; ++i)
        {
        if (n.child[i] >= 0)
          {
          unite_four_aabbs(n.bbox, nodes[n.child[i]].bbox, i);
          }
        else
          {
          n.bbox[0][i] = bb_child[i].bbox_min[0];
          n.bbox[1][i] = bb_child[i].bbox_min[1];
          n.bbox[2][i] = bb_child[i].bbox_min[2];
          n.bbox[3][i] = bb_child[i].bbox_max[0];
          n.bbox[4][i] = bb_child[i].bbox_max[1];
          n.bbox[5][i] = bb_child[i].bbox_max[2];
          }
        }
      int32_t node_id = (int32_t)nodes.size();
      nodes.push_back(n);
      return node_id;
      }
    }
    

  JTKQBVHINLINE qbvh_two_level_with_transformations::qbvh_two_level_with_transformations(const qbvh** objects, const float4x4* transformations, uint32_t nr_of_objects)
    {
    nodes.reserve(nr_of_objects);
    qbvh_voxel total_bb, centroid_bb;
    qbvh_voxel* voxels = new qbvh_voxel[nr_of_objects];
    total_bb.bbox_min = std::numeric_limits<float>::max();
    total_bb.bbox_max = -std::numeric_limits<float>::max();
    centroid_bb.bbox_min = total_bb.bbox_min;
    centroid_bb.bbox_max = total_bb.bbox_max;
    for (uint32_t t = 0; t < nr_of_objects; ++t)
      {
      const qbvh* current_bvh = objects[t];
      unite_four_aabbs(voxels[t], current_bvh->nodes[current_bvh->root_id].bbox);
      transform(voxels[t], transformations[t]);
      voxels[t].centroid = (voxels[t].bbox_min + voxels[t].bbox_max)*0.5f;
      total_bb.bbox_min = min(total_bb.bbox_min, voxels[t].bbox_min);
      total_bb.bbox_max = max(total_bb.bbox_max, voxels[t].bbox_max);
      centroid_bb.bbox_min = min(centroid_bb.bbox_min, voxels[t].centroid);
      centroid_bb.bbox_max = max(centroid_bb.bbox_max, voxels[t].centroid);
      }
    object_ids.resize(nr_of_objects);
    std::iota(object_ids.begin(), object_ids.end(), 0);
    uint32_t sz;
    start = object_ids.begin();
    root_id = construct_tree(sz, voxels, total_bb, centroid_bb, object_ids.begin(), object_ids.end());
    if ((int32_t)root_id < 0)
      {
      qbvh_two_level_node n;
      n.axis0 = 0;
      n.axis1 = 0;
      n.axis2 = 0;
      n.bbox[0] = float4(1.f);
      n.bbox[1] = float4(1.f);
      n.bbox[2] = float4(1.f);
      n.bbox[3] = float4(-1.f);
      n.bbox[4] = float4(-1.f);
      n.bbox[5] = float4(-1.f);

      const auto it = object_ids.begin();
      const auto it_end = it + sz;
      get_bbox(n.bbox, voxels, it, it_end, 0);

      n.child[0] = root_id;
      n.child[1] = std::numeric_limits<int32_t>::min();
      n.child[2] = std::numeric_limits<int32_t>::min();
      n.child[3] = std::numeric_limits<int32_t>::min();
      root_id = 0;
      nodes.push_back(n);
      }
    delete[] voxels;
    }

  JTKQBVHINLINE hit qbvh_two_level_with_transformations::find_closest_triangle(uint32_t& triangle_id, uint32_t& object_id, ray r, const qbvh** objects, const float4x4* inverted_transformations, const vec3<uint32_t>** triangles, const vec3<float>** vertices) const
    {
    hit h;
    h.found = 0;
    h.distance = std::numeric_limits<float>::max();

    vec3<float4> ray_origin(r.orig[0], r.orig[1], r.orig[2]);
    vec3<float4> ray_dir(r.dir[0], r.dir[1], r.dir[2]);

    int32_t ray_dir_sign[3];
    float4 ray_inverse_dir_tmp = reciprocal(r.dir);
    if (fabs(ray_inverse_dir_tmp[0]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[0] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[1]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[1] = std::numeric_limits<float>::max();
    if (fabs(ray_inverse_dir_tmp[2]) == std::numeric_limits<float>::infinity())
      ray_inverse_dir_tmp[2] = std::numeric_limits<float>::max();
    vec3<float4> ray_inverse_dir(ray_inverse_dir_tmp[0], ray_inverse_dir_tmp[1], ray_inverse_dir_tmp[2]);
    float4 t_near(r.t_near);
    float4 t_far(r.t_far);
    ray_dir_sign[0] = r.dir[0] < 0 ? 1 : 0;
    ray_dir_sign[1] = r.dir[1] < 0 ? 1 : 0;
    ray_dir_sign[2] = r.dir[2] < 0 ? 1 : 0;

    std::vector<int32_t> stack;
    stack.reserve(32);
    stack.push_back(root_id);
    while (!stack.empty())
      {
      const auto node = nodes[stack.back()];
      stack.pop_back();
      const uint32_t order_pos = (ray_dir_sign[node.axis0] << 2) | (ray_dir_sign[node.axis1] << 1) | (ray_dir_sign[node.axis2]);
      const uint32_t* local_order = order[order_pos];

      auto isct = intersect(node.bbox, ray_origin, t_near, t_far, ray_dir_sign, ray_inverse_dir);
      if (any(isct))
        {
        const int4 is_leaf = node.child & sign_mask;
        const int contains_at_least_one_leaf = any(is_leaf);
        if (contains_at_least_one_leaf)
          {
          const int4 index = node.child & bit_mask;
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              {
              if (is_leaf[o])
                {
                auto object_it = object_ids.begin() + index[o];
                const uint32_t local_object_id = *object_it;
                uint32_t triangle_candidate;
                ray transformed_r;
                transformed_r.t_far = t_far[o];
                transformed_r.t_near = t_near[o];
                transformed_r.dir = matrix_vector_multiply(inverted_transformations[local_object_id], r.dir);
                transformed_r.orig = matrix_vector_multiply(inverted_transformations[local_object_id], r.orig);
                auto local_hit = objects[local_object_id]->find_closest_triangle(triangle_candidate, transformed_r, triangles[local_object_id], vertices[local_object_id]);
                if (local_hit.found && std::abs(local_hit.distance) < std::abs(h.distance))
                  {
                  h = local_hit;
                  triangle_id = triangle_candidate;
                  object_id = local_object_id;
                  (h.distance > 0) ? t_far = float4(h.distance) : t_near = float4(h.distance);
                  }
                }
              else
                stack.push_back(node.child[o]);
              }
            }
          }
        else
          {
          for (int i = 0; i < 4; ++i)
            {
            const uint32_t o = local_order[i];
            if (isct[o])
              stack.push_back(node.child[o]);
            }
          }
        }
      }

    return h;
    }


  JTKQBVHINLINE int32_t qbvh_two_level_with_transformations::construct_tree(uint32_t& sz, const qbvh_voxel* voxels, qbvh_voxel total_bb, qbvh_voxel centroid_bb, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    assert(first <= last);
    if (first == last)
      {
      sz = 0;
      return std::numeric_limits<int32_t>::min();
      }
    sz = (uint32_t)std::distance(first, last);
    if (sz <= leaf_size)
      {
      uint32_t index = (uint32_t)std::distance(start, first);
      index |= sign_bit;
      return (int32_t)index;
      }
    else
      {
      uint8_t dim;
      qbvh_two_level_node n;
      std::vector<uint32_t>::iterator mid0, mid1, mid2;

      qbvh_voxel bb_left, bb_right, centroid_left, centroid_right;

      qbvh_voxel bb_child[4], centroid_child[4];

      sah_optimized<4>(mid0, bb_left, bb_right, centroid_left, centroid_right, dim, total_bb, centroid_bb, voxels, first, last);
      n.axis0 = dim;
      sah_optimized<4>(mid1, bb_child[0], bb_child[1], centroid_child[0], centroid_child[1], dim, bb_left, centroid_left, voxels, first, mid0);
      n.axis1 = dim;
      sah_optimized<4>(mid2, bb_child[2], bb_child[3], centroid_child[2], centroid_child[3], dim, bb_right, centroid_right, voxels, mid0, last);
      n.axis2 = dim;

      uint32_t sizes[4];

      n.child[0] = construct_tree(sizes[0], voxels, bb_child[0], centroid_child[0], first, mid1);
      n.child[1] = construct_tree(sizes[1], voxels, bb_child[1], centroid_child[1], mid1, mid0);
      n.child[2] = construct_tree(sizes[2], voxels, bb_child[2], centroid_child[2], mid0, mid2);
      n.child[3] = construct_tree(sizes[3], voxels, bb_child[3], centroid_child[3], mid2, last);

      for (int i = 0; i < 4; ++i)
        {
        if (n.child[i] >= 0)
          {
          unite_four_aabbs(n.bbox, nodes[n.child[i]].bbox, i);
          }
        else
          {
          n.bbox[0][i] = bb_child[i].bbox_min[0];
          n.bbox[1][i] = bb_child[i].bbox_min[1];
          n.bbox[2][i] = bb_child[i].bbox_min[2];
          n.bbox[3][i] = bb_child[i].bbox_max[0];
          n.bbox[4][i] = bb_child[i].bbox_max[1];
          n.bbox[5][i] = bb_child[i].bbox_max[2];
          }
        }
      int32_t node_id = (int32_t)nodes.size();
      nodes.push_back(n);
      return node_id;
      }
    }
    
  } // namespace jtk

#endif //#ifndef JTK_QBVH_H


#ifdef JTK_QBVH_IMPLEMENTATION

namespace jtk
  {

  /////////////////////////////////////////////////////////////////////////
  // implementations
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // aligned allocators
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF void* aligned_malloc(size_t size, size_t align)
    {
    if (size == 0)
      return nullptr;

    assert((align & (align - 1)) == 0);
    void* ptr = _mm_malloc(size, align);

    if (ptr == nullptr)
      throw std::bad_alloc();

    return ptr;
    }

  JTKQBVHDEF void aligned_free(void* ptr)
    {
    if (ptr)
      _mm_free(ptr);
    }

  /////////////////////////////////////////////////////////////////////////
  // globals
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // bool4
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF bool all(const bool4& b)
    {
    return _mm_movemask_ps(b.m128) == 0xf;
    }

  JTKQBVHDEF bool any(const bool4& b)
    {
    return _mm_movemask_ps(b.m128) != 0x0;
    }

  JTKQBVHDEF bool none(const bool4& b)
    {
    return _mm_movemask_ps(b.m128) == 0x0;
    }

  JTKQBVHDEF bool4 operator ! (const bool4& a)
    {
    return _mm_xor_ps(a.m128, bool4(true).m128);
    }

  JTKQBVHDEF bool4 operator & (const bool4& a, const bool4& b)
    {
    return _mm_and_ps(a.m128, b.m128);
    }

  JTKQBVHDEF bool4 operator | (const bool4& a, const bool4& b)
    {
    return _mm_or_ps(a.m128, b.m128);
    }

  JTKQBVHDEF bool4 operator ^ (const bool4& a, const bool4& b)
    {
    return _mm_xor_ps(a.m128, b.m128);
    }

  JTKQBVHDEF bool4& operator &= (bool4& a, const bool4& b)
    {
    //return a = a & b; 
    /*
    DANGEROUS: Don't do return a = a & b;
    Lanes will be mixed by compiler in release
    */
    const __m128 res = _mm_and_ps(a.m128, b.m128);
    a.m128 = res;
    return a;
    }

  JTKQBVHDEF bool4& operator |= (bool4& a, const bool4& b)
    {
    //return a = a | b; 
    /*
    DANGEROUS: Don't do return a = a & b;
    Lanes will be mixed by compiler in release
    */
    const __m128 res = _mm_or_ps(a.m128, b.m128);
    a.m128 = res;
    return a;
    }

  JTKQBVHDEF bool4& operator ^= (bool4& a, const bool4& b)
    {
    //return a = a ^ b; 
    /*
    DANGEROUS: Don't do return a = a & b;
    Lanes will be mixed by compiler in release
    */
    const __m128 res = _mm_xor_ps(a.m128, b.m128);
    a.m128 = res;
    return a;
    }

  JTKQBVHDEF bool4 operator != (const bool4& a, const bool4& b)
    {
    return _mm_xor_ps(a.m128, b.m128);
    }

  JTKQBVHDEF bool4 operator == (const bool4& a, const bool4& b)
    {
    return _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_castps_si128(a.m128), _mm_castps_si128(b.m128)));
    }

  /////////////////////////////////////////////////////////////////////////
  // struct int4
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF int4 operator + (const int4& a)
    {
    return a;
    }

  JTKQBVHDEF int4 operator - (const int4& a)
    {
    return _mm_sub_epi32(_mm_set1_epi32(0), a.m128i);
    }

  JTKQBVHDEF int4 operator + (const int4& left, const int4& right)
    {
    return _mm_add_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF int4 operator - (const int4& left, const int4& right)
    {
    return _mm_sub_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF int4 operator * (const int4& left, const int4& right)
    {
    return _mm_mullo_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF int4 operator * (const int4& left, int32_t right)
    {
    return left * int4(right);
    }

  JTKQBVHDEF int4 operator * (int32_t left, const int4& right)
    {
    return int4(left)*right;
    }

  JTKQBVHDEF int4 min(const int4& left, const int4& right)
    {
    return _mm_min_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF int4 max(const int4& left, const int4& right)
    {
    return _mm_max_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF bool4 operator == (const int4& left, const int4& right)
    {
    return _mm_cmpeq_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF bool4 operator != (const int4& left, const int4& right)
    {
    return _mm_andnot_si128(_mm_cmpeq_epi32(left.m128i, right.m128i), not_zero);
    }

  JTKQBVHDEF bool4 operator < (const int4& left, const int4& right)
    {
    return _mm_cmplt_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF bool4 operator > (const int4& left, const int4& right)
    {
    return _mm_cmpgt_epi32(left.m128i, right.m128i);
    }

  JTKQBVHDEF bool4 operator <= (const int4& left, const int4& right)
    {
    return _mm_andnot_si128(_mm_cmpgt_epi32(left.m128i, right.m128i), not_zero);
    }

  JTKQBVHDEF bool4 operator >= (const int4& left, const int4& right)
    {
    return _mm_andnot_si128(_mm_cmplt_epi32(left.m128i, right.m128i), not_zero);
    }

  JTKQBVHDEF int4 masked_update(const bool4& mask, const int4& original, const int4& updated_values)
    {
    return _mm_castps_si128(_mm_or_ps(_mm_and_ps(mask.m128, _mm_castsi128_ps(updated_values.m128i)), _mm_andnot_ps(mask.m128, _mm_castsi128_ps(original.m128i))));
    }

  JTKQBVHDEF int4 operator & (const int4& left, const int4& right)
    {
    return _mm_and_si128(left.m128i, right.m128i);
    }

  JTKQBVHDEF int4 operator | (const int4& left, const int4& right)
    {
    return _mm_and_si128(left.m128i, right.m128i);
    }

#ifndef _JTK_FOR_ARM
  JTKQBVHDEF int4 operator >> (const int4& a, int n)
    {
    return _mm_srai_epi32(a.m128i, n);
    }

  JTKQBVHDEF int4 operator << (const int4& a, int n)
    {
    return _mm_slli_epi32(a.m128i, n);
    }
#endif

  JTKQBVHDEF int any(const int4& a)
    {
    return _mm_movemask_ps(_mm_castsi128_ps(a.m128i));
    }

  /////////////////////////////////////////////////////////////////////////
  // struct float4
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF float4 operator + (const float4& a)
    {
    return a;
    }

  JTKQBVHDEF float4 operator - (const float4& a)
    {
    return _mm_xor_ps(a.m128, _mm_castsi128_ps(_mm_set1_epi32(0x80000000)));
    }

  JTKQBVHDEF float4 operator + (const float4& left, const float4& right)
    {
    return _mm_add_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 operator - (const float4& left, const float4& right)
    {
    return _mm_sub_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 operator * (const float4& left, const float4& right)
    {
    return _mm_mul_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 operator * (const float4& left, float right)
    {
    return left * float4(right);
    }

  JTKQBVHDEF float4 operator * (float left, const float4& right)
    {
    return float4(left)*right;
    }

  JTKQBVHDEF float4 operator / (const float4& left, const float4& right)
    {
    return _mm_div_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 operator / (const float4& left, float right)
    {
    return left / float4(right);
    }

  JTKQBVHDEF float4 operator / (float left, const float4& right)
    {
    return float4(left) / right;
    }

  JTKQBVHDEF float4 min(const float4& left, const float4& right)
    {
    return _mm_min_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 max(const float4& left, const float4& right)
    {
    return _mm_max_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float min_horizontal(const float4& x)
    {
    __m128 max1 = _mm_shuffle_ps(x.m128, x.m128, _MM_SHUFFLE(0, 0, 3, 2));
    __m128 max2 = _mm_min_ps(x.m128, max1);
    __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128 max4 = _mm_min_ps(max2, max3);
    float result = _mm_cvtss_f32(max4);
    return result;
    }

  JTKQBVHDEF float max_horizontal(const float4& x)
    {
    __m128 max1 = _mm_shuffle_ps(x.m128, x.m128, _MM_SHUFFLE(0, 0, 3, 2));
    __m128 max2 = _mm_max_ps(x.m128, max1);
    __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128 max4 = _mm_max_ps(max2, max3);
    float result = _mm_cvtss_f32(max4);
    return result;
    }

  JTKQBVHDEF float4 cross(const float4& left, const float4& right)
    {
    float4 rs(_mm_shuffle_ps(right.m128, right.m128, _MM_SHUFFLE(3, 0, 2, 1)));
    float4 ls(_mm_shuffle_ps(left.m128, left.m128, _MM_SHUFFLE(3, 0, 2, 1)));
    float4 res = left * rs - ls * right;
    return float4(_mm_shuffle_ps(res.m128, res.m128, _MM_SHUFFLE(3, 0, 2, 1)));
    }

  JTKQBVHDEF float dot(const float4& left, const float4& right)
    {
    return _mm_cvtss_f32(_mm_dp_ps(left.m128, right.m128, 0x7F));
    }

  JTKQBVHDEF float dot4(const float4& left, const float4& right)
    {
    return _mm_cvtss_f32(_mm_dp_ps(left.m128, right.m128, 255));
    }

  JTKQBVHDEF float4 abs(const float4& a)
    {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return _mm_and_ps(a.m128, mask);
    }

  JTKQBVHDEF float4 sqrt(const float4& a)
    {
    return _mm_sqrt_ps(a.m128);
    }

  JTKQBVHDEF float4 rsqrt(const float4& a)
    {
    __m128 mask = _mm_cmpeq_ps(_mm_set1_ps(0.f), a.m128);
    __m128 res = _mm_rsqrt_ps(a.m128);
    __m128 muls = _mm_mul_ps(_mm_mul_ps(a.m128, res), res);
    auto res_newton_raphson = _mm_mul_ps(_mm_mul_ps(half_128, res), _mm_sub_ps(three_128, muls));
    return _mm_or_ps(_mm_and_ps(mask, res), _mm_andnot_ps(mask, res_newton_raphson));
    }

  JTKQBVHDEF float4 reciprocal(const float4& a)
    {
    __m128 mask = _mm_cmpeq_ps(_mm_set1_ps(0.f), a.m128);
    auto res = _mm_rcp_ps(a.m128);
    __m128 muls = _mm_mul_ps(a.m128, _mm_mul_ps(res, res));
    auto res_newton_raphson = _mm_sub_ps(_mm_add_ps(res, res), muls);
    return _mm_or_ps(_mm_and_ps(mask, res), _mm_andnot_ps(mask, res_newton_raphson));
    }

  JTKQBVHDEF bool4 operator == (const float4& left, const float4& right)
    {
    return _mm_cmpeq_ps(left.m128, right.m128);
    }

  JTKQBVHDEF bool4 operator != (const float4& left, const float4& right)
    {
    return _mm_cmpneq_ps(left.m128, right.m128);
    }

  JTKQBVHDEF bool4 operator < (const float4& left, const float4& right)
    {
    return _mm_cmplt_ps(left.m128, right.m128);
    }

  JTKQBVHDEF bool4 operator > (const float4& left, const float4& right)
    {
    return _mm_cmpnle_ps(left.m128, right.m128);
    }

  JTKQBVHDEF bool4 operator <= (const float4& left, const float4& right)
    {
    return _mm_cmple_ps(left.m128, right.m128);
    }

  JTKQBVHDEF bool4 operator >= (const float4& left, const float4& right)
    {
    return _mm_cmpnlt_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 unpacklo(const float4& left, const float4& right)
    {
    return _mm_unpacklo_ps(left.m128, right.m128);
    }

  JTKQBVHDEF float4 unpackhi(const float4& left, const float4& right)
    {
    return _mm_unpackhi_ps(left.m128, right.m128);
    }

  JTKQBVHDEF void transpose(float4& r0, float4& r1, float4& r2, float4& r3, const float4& c0, const float4& c1, const float4& c2, const float4& c3)
    {
    float4 l02(unpacklo(c0.m128, c2.m128));
    float4 h02(unpackhi(c0.m128, c2.m128));
    float4 l13(unpacklo(c1.m128, c3.m128));
    float4 h13(unpackhi(c1.m128, c3.m128));
    r0 = unpacklo(l02, l13);
    r1 = unpackhi(l02, l13);
    r2 = unpacklo(h02, h13);
    r3 = unpackhi(h02, h13);
    }

  JTKQBVHDEF float4 masked_update(const bool4& mask, const float4& original, const float4& updated_values)
    {
    return _mm_or_ps(_mm_and_ps(mask.m128, updated_values.m128), _mm_andnot_ps(mask.m128, original.m128));
    }

  JTKQBVHDEF float4 masked_update(const int4& mask, const float4& original, const float4& updated_values)
    {
    const __m128 m = _mm_castsi128_ps(mask.m128i);
    return _mm_or_ps(_mm_and_ps(m, updated_values.m128), _mm_andnot_ps(m, original.m128));
    }

  /////////////////////////////////////////////////////////////////////////
  // struct float4x4
  /////////////////////////////////////////////////////////////////////////

  // COLUMN MAJOR 4x4 MATRIX

  JTKQBVHDEF float4x4 get_identity()
    {
    float4x4 m(_mm_set_ps(0.f, 0.f, 0.f, 1.f), _mm_set_ps(0.f, 0.f, 1.f, 0.f), _mm_set_ps(0.f, 1.f, 0.f, 0.f), _mm_set_ps(1.f, 0.f, 0.f, 0.f));
    return m;
    }

  JTKQBVHDEF float4x4 make_translation(float x, float y, float z)
    {
    float4x4 m(_mm_set_ps(0.f, 0.f, 0.f, 1.f), _mm_set_ps(0.f, 0.f, 1.f, 0.f), _mm_set_ps(0.f, 1.f, 0.f, 0.f), _mm_set_ps(1.f, z, y, x));
    return m;
    }

  JTKQBVHDEF void solve_roll_pitch_yaw_transformation(float& rx, float& ry, float& rz, float& tx, float& ty, float& tz, const float4x4& m)
    {
    rz = std::atan2(m.col[0][1], m.col[0][0]);
    const auto sg = std::sin(rz);
    const auto cg = std::cos(rz);    
    ry = std::atan2(-m.col[0][2], m.col[0][0] * cg + m.col[0][1] * sg);    
    rx = std::atan2(m.col[2][0] * sg - m.col[2][1] * cg, m.col[1][1] * cg - m.col[1][0] * sg);    
    tx = m.col[3][0];    
    ty = m.col[3][1];    
    tz = m.col[3][2];
    }

  JTKQBVHDEF float4x4 compute_from_roll_pitch_yaw_transformation(
    float rx, float ry, float rz,
    float tx, float ty, float tz)
    {
    float4x4 m = get_identity();
    float ca = std::cos(rx);
    float sa = std::sin(rx);
    float cb = std::cos(ry);
    float sb = std::sin(ry);
    float cg = std::cos(rz);
    float sg = std::sin(rz);    
    m.col[0][0] = cb * cg;
    m.col[1][0] = cg * sa * sb - ca * sg;
    m.col[2][0] = sa * sg + ca * cg * sb;
    m.col[0][1] = cb * sg;
    m.col[1][1] = sa * sb * sg + ca * cg;
    m.col[2][1] = ca * sb * sg - cg * sa;
    m.col[0][2] = -sb;
    m.col[1][2] = cb * sa;
    m.col[2][2] = ca * cb;
    m.col[3][0] = tx;
    m.col[3][1] = ty;
    m.col[3][2] = tz;
    return m;
    }

  JTKQBVHDEF float4x4 quaternion_to_rotation(const float4& quaternion)
    {
    jtk::float4x4 rot;
    rot[0] = 1.f - 2.f * (quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2]);
    rot[1] = 2.f * (quaternion[0] * quaternion[1] - quaternion[2] * quaternion[3]);
    rot[2] = 2.f * (quaternion[2] * quaternion[0] + quaternion[1] * quaternion[3]);
    rot[3] = 0.f;

    rot[4] = 2.f * (quaternion[0] * quaternion[1] + quaternion[2] * quaternion[3]);
    rot[5] = 1.f - 2.f * (quaternion[2] * quaternion[2] + quaternion[0] * quaternion[0]);
    rot[6] = 2.f * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3]);
    rot[7] = 0.f;

    rot[8] = 2.f * (quaternion[2] * quaternion[0] - quaternion[1] * quaternion[3]);
    rot[9] = 2.f * (quaternion[1] * quaternion[2] + quaternion[0] * quaternion[3]);
    rot[10] = 1.f - 2.f * (quaternion[1] * quaternion[1] + quaternion[0] * quaternion[0]);
    rot[11] = 0.f;

    rot[12] = 0.f;
    rot[13] = 0.f;
    rot[14] = 0.f;
    rot[15] = 1.f;
    return rot;
    }

  JTKQBVHDEF float4 roll_pitch_yaw_to_quaternion(float rx, float ry, float rz)
    {    
    float cy = std::cos(rz * 0.5f);
    float sy = std::sin(rz * 0.5f);
    float cp = std::cos(ry * 0.5f);
    float sp = std::sin(ry * 0.5f);
    float cr = std::cos(rx * 0.5f);
    float sr = std::sin(rx * 0.5f);

    float4 q;    
    q[0] = sr * cp * cy - cr * sp * sy;
    q[1] = cr * sp * cy + sr * cp * sy;
    q[2] = cr * cp * sy - sr * sp * cy;
    q[3] = cr * cp * cy + sr * sp * sy;
    return q;
    }

  JTKQBVHDEF float4x4 transpose(const float4x4& m)
    {
    float4x4 out;
    transpose(out.col[0], out.col[1], out.col[2], out.col[3], m.col[0], m.col[1], m.col[2], m.col[3]);
    return out;
    }

  JTKQBVHDEF float4x4 invert_orthonormal(const float4x4& m)
    {
    float4x4 out;
    transpose(out.col[0], out.col[1], out.col[2], out.col[3], m.col[0], m.col[1], m.col[2], _mm_set_ps(1.f, 0.f, 0.f, 0.f));
    out.col[3] = -(out.col[0] * m[12] + out.col[1] * m[13] + out.col[2] * m[14]);
    out.f[15] = 1.f;
    return out;
    }

  // for column major matrix
  // we use __m128 to represent 2x2 matrix as A = | A0  A1 |
  //                                              | A2  A3 |
  // 2x2 column major matrix multiply A*B
  JTKQBVHDEF __m128 mat2mul(__m128 vec1, __m128 vec2)
    {
    const auto vec3 = _mm_mul_ps(vec1, _mm_shuffle_ps(vec2, vec2, _MM_SHUFFLE(3, 3, 0, 0)));
    const auto vec4 = _mm_mul_ps(_mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(1, 0, 3, 2)), _mm_shuffle_ps(vec2, vec2, _MM_SHUFFLE(2, 2, 1, 1)));
    return _mm_add_ps(vec3, vec4);
    }

  // 2x2 column major matrix adjugate multiply (A#)*B
  JTKQBVHDEF __m128 mat2adjmul(__m128 vec1, __m128 vec2)
    {
    const auto vec3 = _mm_mul_ps(_mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(0, 3, 0, 3)), vec2);
    const auto vec4 = _mm_mul_ps(_mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(1, 2, 1, 2)), _mm_shuffle_ps(vec2, vec2, _MM_SHUFFLE(2, 3, 0, 1)));
    return _mm_sub_ps(vec3, vec4);
    }

  // 2x2 column major matrix multiply adjugate A*(B#)
  JTKQBVHDEF __m128 mat2muladj(__m128 vec1, __m128 vec2)
    {
    const auto vec3 = _mm_mul_ps(vec1, _mm_shuffle_ps(vec2, vec2, _MM_SHUFFLE(0, 0, 3, 3)));
    const auto vec4 = _mm_mul_ps(_mm_shuffle_ps(vec1, vec1, _MM_SHUFFLE(1, 0, 3, 2)), _mm_shuffle_ps(vec2, vec2, _MM_SHUFFLE(2, 2, 1, 1)));
    return _mm_sub_ps(vec3, vec4);
    }

  JTKQBVHDEF float4x4 invert(const float4x4& m)
    {
    float4x4 out;
    // sub matrices
    __m128 A = _mm_shuffle_ps(m.col[0].m128, m.col[1].m128, _MM_SHUFFLE(1, 0, 1, 0));
    __m128 C = _mm_shuffle_ps(m.col[0].m128, m.col[1].m128, _MM_SHUFFLE(3, 2, 3, 2));
    __m128 B = _mm_shuffle_ps(m.col[2].m128, m.col[3].m128, _MM_SHUFFLE(1, 0, 1, 0));
    __m128 D = _mm_shuffle_ps(m.col[2].m128, m.col[3].m128, _MM_SHUFFLE(3, 2, 3, 2));

    // determinant as (|A| |B| |C| |D|)
    __m128 detSub = _mm_sub_ps(_mm_mul_ps(
      _mm_shuffle_ps(m.col[0].m128, m.col[2].m128, _MM_SHUFFLE(2, 0, 2, 0)),
      _mm_shuffle_ps(m.col[1].m128, m.col[3].m128, _MM_SHUFFLE(3, 1, 3, 1))),
      _mm_mul_ps(_mm_shuffle_ps(m.col[0].m128, m.col[2].m128, _MM_SHUFFLE(3, 1, 3, 1)),
        _mm_shuffle_ps(m.col[1].m128, m.col[3].m128, _MM_SHUFFLE(2, 0, 2, 0))));
    __m128 detA = _mm_shuffle_ps(detSub, detSub, _MM_SHUFFLE(0, 0, 0, 0));
    __m128 detC = _mm_shuffle_ps(detSub, detSub, _MM_SHUFFLE(1, 1, 1, 1));
    __m128 detB = _mm_shuffle_ps(detSub, detSub, _MM_SHUFFLE(2, 2, 2, 2));
    __m128 detD = _mm_shuffle_ps(detSub, detSub, _MM_SHUFFLE(3, 3, 3, 3));

    // let iM = 1/|M| * | X  Y |
    //                  | Z  W |

    // D#C
    __m128 D_C = mat2adjmul(D, C);
    // A#B
    __m128 A_B = mat2adjmul(A, B);
    // X# = |D|A - B(D#C)
    __m128 X_ = _mm_sub_ps(_mm_mul_ps(detD, A), mat2mul(B, D_C));
    // W# = |A|D - C(A#B)
    __m128 W_ = _mm_sub_ps(_mm_mul_ps(detA, D), mat2mul(C, A_B));

    // |M| = |A|*|D| + ... (continue later)
    __m128 detM = _mm_mul_ps(detA, detD);

    // Y# = |B|C - D(A#B)#
    __m128 Y_ = _mm_sub_ps(_mm_mul_ps(detB, C), mat2muladj(D, A_B));
    // Z# = |C|B - A(D#C)#
    __m128 Z_ = _mm_sub_ps(_mm_mul_ps(detC, B), mat2muladj(A, D_C));

    // |M| = |A|*|D| + |B|*|C| ... (continue later)
    detM = _mm_add_ps(detM, _mm_mul_ps(detB, detC));

    // tr((A#B)(D#C))
    __m128 tr = _mm_mul_ps(A_B, _mm_shuffle_ps(D_C, D_C, _MM_SHUFFLE(3, 1, 2, 0)));
    tr = _mm_hadd_ps(tr, tr);
    tr = _mm_hadd_ps(tr, tr);
    // |M| = |A|*|D| + |B|*|C| - tr((A#B)(D#C)
    detM = _mm_sub_ps(detM, tr);

    const __m128 adjSignMask = _mm_setr_ps(1.f, -1.f, -1.f, 1.f);
    // (1/|M|, -1/|M|, -1/|M|, 1/|M|)
    __m128 rDetM = _mm_div_ps(adjSignMask, detM);

    X_ = _mm_mul_ps(X_, rDetM);
    Y_ = _mm_mul_ps(Y_, rDetM);
    Z_ = _mm_mul_ps(Z_, rDetM);
    W_ = _mm_mul_ps(W_, rDetM);

    // apply adjugate and store, here we combine adjugate shuffle and store shuffle
    out.col[0] = float4(_mm_shuffle_ps(X_, Z_, _MM_SHUFFLE(1, 3, 1, 3)));
    out.col[1] = float4(_mm_shuffle_ps(X_, Z_, _MM_SHUFFLE(0, 2, 0, 2)));
    out.col[2] = float4(_mm_shuffle_ps(Y_, W_, _MM_SHUFFLE(1, 3, 1, 3)));
    out.col[3] = float4(_mm_shuffle_ps(Y_, W_, _MM_SHUFFLE(0, 2, 0, 2)));

    return out;
    }

  JTKQBVHDEF float4 matrix_vector_multiply(const float4x4& m, const float4& v)
    {
    float4 out = m.col[0] * v[0] + m.col[1] * v[1] + m.col[2] * v[2] + m.col[3] * v[3];
    return out;
    }

  JTKQBVHDEF float4x4 matrix_matrix_multiply(const float4x4& left, const float4x4& right)
    {
    float4x4 out;
    float4 r[4];
    transpose(r[0], r[1], r[2], r[3], left.col[0], left.col[1], left.col[2], left.col[3]);
    out[0] = _mm_cvtss_f32(_mm_dp_ps(r[0].m128, right.col[0].m128, 255));
    out[1] = _mm_cvtss_f32(_mm_dp_ps(r[1].m128, right.col[0].m128, 255));
    out[2] = _mm_cvtss_f32(_mm_dp_ps(r[2].m128, right.col[0].m128, 255));
    out[3] = _mm_cvtss_f32(_mm_dp_ps(r[3].m128, right.col[0].m128, 255));
    out[4] = _mm_cvtss_f32(_mm_dp_ps(r[0].m128, right.col[1].m128, 255));
    out[5] = _mm_cvtss_f32(_mm_dp_ps(r[1].m128, right.col[1].m128, 255));
    out[6] = _mm_cvtss_f32(_mm_dp_ps(r[2].m128, right.col[1].m128, 255));
    out[7] = _mm_cvtss_f32(_mm_dp_ps(r[3].m128, right.col[1].m128, 255));
    out[8] = _mm_cvtss_f32(_mm_dp_ps(r[0].m128, right.col[2].m128, 255));
    out[9] = _mm_cvtss_f32(_mm_dp_ps(r[1].m128, right.col[2].m128, 255));
    out[10] = _mm_cvtss_f32(_mm_dp_ps(r[2].m128, right.col[2].m128, 255));
    out[11] = _mm_cvtss_f32(_mm_dp_ps(r[3].m128, right.col[2].m128, 255));
    out[12] = _mm_cvtss_f32(_mm_dp_ps(r[0].m128, right.col[3].m128, 255));
    out[13] = _mm_cvtss_f32(_mm_dp_ps(r[1].m128, right.col[3].m128, 255));
    out[14] = _mm_cvtss_f32(_mm_dp_ps(r[2].m128, right.col[3].m128, 255));
    out[15] = _mm_cvtss_f32(_mm_dp_ps(r[3].m128, right.col[3].m128, 255));
    return out;
    }

  JTKQBVHDEF float4x4 operator + (const float4x4& left, const float4x4& right)
    {
    return float4x4(left.col[0] + right.col[0], left.col[1] + right.col[1], left.col[2] + right.col[2], left.col[3] + right.col[3]);
    }

  JTKQBVHDEF float4x4 operator - (const float4x4& left, const float4x4& right)
    {
    return float4x4(left.col[0] - right.col[0], left.col[1] - right.col[1], left.col[2] - right.col[2], left.col[3] - right.col[3]);
    }

  JTKQBVHDEF float4x4 operator / (const float4x4& left, float value)
    {
    return float4x4(left.col[0] / value, left.col[1] / value, left.col[2] / value, left.col[3] / value);
    }

  JTKQBVHDEF float4x4 operator * (const float4x4& left, float value)
    {
    return float4x4(left.col[0] * value, left.col[1] * value, left.col[2] * value, left.col[3] * value);
    }

  JTKQBVHDEF float4x4 operator * (float value, const float4x4& right)
    {
    return float4x4(right.col[0] * value, right.col[1] * value, right.col[2] * value, right.col[3] * value);
    }


  /////////////////////////////////////////////////////////////////////////
  // intersection
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF woop_precompute intersect_woop_precompute(const float4& r_dir)
    {
    woop_precompute out;
    const auto abs_dir = abs(r_dir);
    out.kz = 2;
    if (abs_dir[0] > abs_dir[1])
      {
      if (abs_dir[0] > abs_dir[2])
        out.kz = 0;
      }
    else
      {
      if (abs_dir[1] > abs_dir[2])
        out.kz = 1;
      }
    out.kx = modulo[out.kz + 1];
    out.ky = modulo[out.kx + 1];

    if (r_dir[out.kz] < (typename float4::value_type)0)
      {
      const uint32_t tmp = out.kx;
      out.kx = out.ky;
      out.ky = tmp;
      }

    out.Sz = (typename float4::value_type)1 / r_dir[out.kz];
    out.Sx = r_dir[out.kx] * out.Sz;
    out.Sy = r_dir[out.ky] * out.Sz;

    return out;
    }

  JTKQBVHDEF void intersect_woop(const woop_triangle& acc, const woop_precompute& pre, const vec3<float4>& r_orig, const float4& t_near, const float4& t_far, hit4& h)
    {
    h.found = bool4(true);

    const auto A = acc.v0 - r_orig;
    const auto B = acc.v1 - r_orig;
    const auto C = acc.v2 - r_orig;

    const float4 Ax = A[pre.kx] - pre.Sx*A[pre.kz];
    const float4 Ay = A[pre.ky] - pre.Sy*A[pre.kz];
    const float4 Bx = B[pre.kx] - pre.Sx*B[pre.kz];
    const float4 By = B[pre.ky] - pre.Sy*B[pre.kz];
    const float4 Cx = C[pre.kx] - pre.Sx*C[pre.kz];
    const float4 Cy = C[pre.ky] - pre.Sy*C[pre.kz];

    float4 U = Cx * By - Cy * Bx;
    float4 V = Ax * Cy - Ay * Cx;
    float4 W = Bx * Ay - By * Ax;

    h.found &= (((U <= float4(0)) & (V <= float4(0)) & (W <= float4(0))) | ((U >= float4(0)) & (V >= float4(0)) & (W >= float4(0))));

    if (none(h.found))
      return;

    const float4 det = U + V + W;

    h.found &= (det != float4(0));

    if (none(h.found))
      return;

    const float4 inv_det = reciprocal(det);

    const float4 Az = pre.Sz*A[pre.kz];
    const float4 Bz = pre.Sz*B[pre.kz];
    const float4 Cz = pre.Sz*C[pre.kz];
    const float4 T = U * Az + V * Bz + W * Cz;
    const float4 t = T * inv_det;

    h.found &= ((t_far > t) & (t > t_near));

    h.u = V * inv_det;
    h.v = W * inv_det;
    h.distance = t;
    }

  JTKQBVHDEF bool4 intersect(const float4* aabb, const vec3<float4>& r_origin, const float4& t_near, const float4& t_far, const int32_t* ray_dir_sign, const vec3<float4>& ray_inverse_dir)
    {
    bool4 res(true);
    const float4 min_x = ray_dir_sign[0] ? aabb[3] : aabb[0];
    const float4 max_x = ray_dir_sign[0] ? aabb[0] : aabb[3];
    const float4 min_y = ray_dir_sign[1] ? aabb[4] : aabb[1];
    const float4 max_y = ray_dir_sign[1] ? aabb[1] : aabb[4];
    const float4 min_z = ray_dir_sign[2] ? aabb[5] : aabb[2];
    const float4 max_z = ray_dir_sign[2] ? aabb[2] : aabb[5];

    float4 txmin = (min_x - r_origin[0]) * ray_inverse_dir[0];
    float4 txmax = (max_x - r_origin[0]) * ray_inverse_dir[0];

    const float4 tymin = (min_y - r_origin[1]) * ray_inverse_dir[1];
    const float4 tymax = (max_y - r_origin[1]) * ray_inverse_dir[1];

    res &= (txmin <= tymax);
    res &= (tymin <= txmax);

    if (none(res))
      return res;

    const float4 tzmin = (min_z - r_origin[2]) * ray_inverse_dir[2];
    const float4 tzmax = (max_z - r_origin[2]) * ray_inverse_dir[2];

    txmin = max(txmin, tymin);
    txmax = min(txmax, tymax);
    res &= (txmin <= tzmax);
    res &= (tzmin <= txmax);

    if (none(res))
      return res;

    txmin = max(txmin, tzmin);
    txmax = min(txmax, tzmax);
    res &= (txmin <= t_far);
    res &= (txmax >= t_near);
    return res;
    }

  JTKQBVHDEF void intersect_sphere(const vec3<float4>& sphere_origin, const float4& sphere_radius, const vec3<float4>& r_orig, const vec3<float4>& ray_dir, const float4& t_near, const float4& t_far, spherehit4& h)
    {
#if 0
    h.found = bool4(true);
    const vec3<float4> L = sphere_origin - r_orig;
    const float4 tca = dot(L, ray_dir);
    const float4 d2 = dot(L, L) - tca * tca;
    const float4 sphere_radius_sqr = sphere_radius * sphere_radius;
    h.found &= (d2 <= sphere_radius_sqr);
    if (none(h.found))
      return;
    const float4 thc = sqrt(sphere_radius_sqr - d2);
    const float4 t0 = tca - thc;
    const float4 t1 = tca + thc;
    const bool4 t0_is_valid = (t_near <= t0) & (t0 <= t_far);
    const bool4 t1_is_valid = (t_near <= t1) & (t1 <= t_far);
    h.found &= t0_is_valid | t1_is_valid;
    if (none(h.found))
      return;
    const bool4 t0_is_closest = abs(t0) < abs(t1);
    const float4 t_closest = masked_update(t0_is_closest, t1, t0);
    h.distance = masked_update(t0_is_valid&t1_is_valid, masked_update(t0_is_valid, t1, t0), t_closest);
#else
    /*
    Precision Improvements for Ray/Sphere Intersection
    Authors:
      Eric Haines (NVIDIA)
      Johannes Gunther (Intel)
      Tomas Akenine-Möller
    Publication Date:
      Friday, March 1, 2019
    Published in:
      Ray Tracing Gems
    */

    /*
    h.found = bool4(true);
    const vec3<float4> f = sphere_origin - r_orig;
    const float4 sphere_radius_sqr = sphere_radius * sphere_radius;
    const vec3<float4> l = f - dot(f, ray_dir)*ray_dir;
    const float4 discr = 4.f*(sphere_radius_sqr - dot(l, l));
    h.found &= discr >= 0.f;
    if (none(h.found))
      return;
    const float4 b = dot(f, ray_dir);
    const float4 sqrt_discr_div_2 = sqrt(discr) / 2.f;
    const float4 t0 = b - sqrt_discr_div_2;
    const float4 t1 = b + sqrt_discr_div_2;
    const bool4 t0_is_valid = (t_near <= t0) & (t0 <= t_far);
    const bool4 t1_is_valid = (t_near <= t1) & (t1 <= t_far);
    h.found &= t0_is_valid | t1_is_valid;
    if (none(h.found))
      return;
    const bool4 t0_is_closest = abs(t0) < abs(t1);
    const float4 t_closest = masked_update(t0_is_closest, t1, t0);
    h.distance = masked_update(t0_is_valid&t1_is_valid, masked_update(t0_is_valid, t1, t0), t_closest);
    */

    h.found = bool4(true);
    const vec3<float4> f = sphere_origin - r_orig;
    const float4 b = dot(f, ray_dir);
    const float4 sphere_radius_sqr = sphere_radius * sphere_radius;
    const vec3<float4> l = b * ray_dir - f;
    const float4 delta = sphere_radius_sqr - dot(l, l);
    h.found &= delta >= 0.f;
    if (none(h.found))
      return;
    const float4 c = dot(f, f) - sphere_radius_sqr;
    const float4 q = b + masked_update(b > float4(0.f), float4(-1.f), float4(1.f))*sqrt(delta);
    const float4 t0 = c / q;
    const float4 t1 = q;
    const bool4 t0_is_valid = (t_near <= t0) & (t0 <= t_far);
    const bool4 t1_is_valid = (t_near <= t1) & (t1 <= t_far);
    h.found &= t0_is_valid | t1_is_valid;
    if (none(h.found))
      return;
    const bool4 t0_is_closest = abs(t0) < abs(t1);
    const float4 t_closest = masked_update(t0_is_closest, t1, t0);
    h.distance = masked_update(t0_is_valid&t1_is_valid, masked_update(t0_is_valid, t1, t0), t_closest);
#endif
    }

  /////////////////////////////////////////////////////////////////////////
  // distance
  /////////////////////////////////////////////////////////////////////////

  //http://jcgt.org/published/0003/04/05/paper.pdf
  JTKQBVHDEF void distance_sqr(const woop_triangle& acc, const vec3<float4>& point, distance4& dist)
    {
    const vec3<float4> ab = acc.v1 - acc.v0;
    const vec3<float4> ac = acc.v2 - acc.v0;
    const vec3<float4> ap = point - acc.v0;
    const float4 d1 = dot(ab, ap);
    const float4 d2 = dot(ac, ap);
    const auto mask1 = (d1 <= float4(0)) & (d2 <= float4(0));
    // closest point is v0
    dist.u = float4(0);
    dist.v = float4(0);
    auto exit(mask1);
    if (all(exit))
      {
      dist.distance_sqr = length_sqr(acc.v0 - point);
      return;
      }

    const vec3<float4> bp = ap - ab;
    const float4 d3 = dot(ab, bp);
    const float4 d4 = dot(ac, bp);
    const auto mask2 = (d3 >= float4(0)) & (d4 <= d3);
    // closest point is v1  
    dist.u = masked_update(exit, masked_update(mask2, dist.u, float4(1)), dist.u);
    dist.v = masked_update(exit, masked_update(mask2, dist.v, float4(0)), dist.v);
    exit |= mask2;
    if (all(exit))
      {
      const vec3<float4> closest_point = acc.v0 + dist.u*ab + dist.v*ac;
      dist.distance_sqr = length_sqr(closest_point - point);
      return;
      }

    const vec3<float4> cp = ap - ac;
    const float4 d5 = dot(ab, cp);
    const float4 d6 = dot(ac, cp);
    const auto mask3 = (d6 >= float4(0)) & (d5 <= d6);
    // closest point is v2  
    dist.u = masked_update(exit, masked_update(mask3, dist.u, float4(0)), dist.u);
    dist.v = masked_update(exit, masked_update(mask3, dist.v, float4(1)), dist.v);
    exit |= mask3;
    if (all(exit))
      {
      const vec3<float4> closest_point = acc.v0 + dist.u*ab + dist.v*ac;
      dist.distance_sqr = length_sqr(closest_point - point);
      return;
      }

    const float4 vc = d1 * d4 - d3 * d2;
    const auto mask4 = (vc <= float4(0)) & (d1 >= float4(0)) & (d3 <= float4(0));
    const float4 v1 = d1 / (d1 - d3);
    //const vec3<float4> answer1 = acc.v0 + v1 * ab;
    // closest point is on the line ab  
    dist.u = masked_update(exit, masked_update(mask4, dist.u, v1), dist.u);
    dist.v = masked_update(exit, masked_update(mask4, dist.v, 0.0), dist.v);
    exit |= mask4;
    if (all(exit))
      {
      const vec3<float4> closest_point = acc.v0 + dist.u*ab + dist.v*ac;
      dist.distance_sqr = length_sqr(closest_point - point);
      return;
      }

    const float4 vb = d5 * d2 - d1 * d6;
    const auto mask5 = (vb <= float4(0)) & (d2 >= float4(0)) & (d6 <= float4(0));
    const float4 w1 = d2 / (d2 - d6);
    //const vec3<float4> answer2 = acc.v0 + w1 * ac;
    // closest point is on the line ac  
    dist.u = masked_update(exit, masked_update(mask5, dist.u, float4(0)), dist.u);
    dist.v = masked_update(exit, masked_update(mask5, dist.v, w1), dist.v);
    exit |= mask5;
    if (all(exit))
      {
      const vec3<float4> closest_point = acc.v0 + dist.u*ab + dist.v*ac;
      dist.distance_sqr = length_sqr(closest_point - point);
      return;
      }

    const float4 va = d3 * d6 - d5 * d4;
    const auto mask6 = (va <= float4(0)) & ((d4 - d3) >= float4(0)) & ((d5 - d6) >= float4(0));
    const float4 w2 = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    //const vec3<float4> answer3 = (float4(1) - w2)*ab + w2 * ac;
    // closest point is on the line bc  
    dist.u = masked_update(exit, masked_update(mask6, dist.u, float4(1) - w2), dist.u);
    dist.v = masked_update(exit, masked_update(mask6, dist.v, w2), dist.v);
    exit |= mask6;
    if (all(exit))
      {
      const vec3<float4> closest_point = acc.v0 + dist.u*ab + dist.v*ac;
      dist.distance_sqr = length_sqr(closest_point - point);
      return;
      }

    const float4 denom = float4(1) / (va + vb + vc);
    const float4 v2 = vb * denom;
    const float4 w3 = vc * denom;
    //const vec3<float4> answer4 = acc.v0 + ab * v2 + ac * w3;
    const auto mask7 = !exit;
    dist.u = masked_update(exit, masked_update(mask7, dist.u, v2), dist.u);
    dist.v = masked_update(exit, masked_update(mask7, dist.v, w3), dist.v);
    const vec3<float4> closest_point = acc.v0 + dist.u*ab + dist.v*ac;
    dist.distance_sqr = length_sqr(closest_point - point);
    }

  JTKQBVHDEF float4 distance_sqr(const float4* aabb, const vec3<float4>& point)
    {
    const float4 x = point[0] - min(aabb[3], max(aabb[0], point[0]));
    const float4 y = point[1] - min(aabb[4], max(aabb[1], point[1]));
    const float4 z = point[2] - min(aabb[5], max(aabb[2], point[2]));

    return x * x + y * y + z * z;
    }

  /////////////////////////////////////////////////////////////////////////
  // transformation
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF float4x4 make_identity()
    {
    return get_identity();
    }

  JTKQBVHDEF vec3<float> transform(const float4x4& matrix, const vec3<float>& pt)
    {
    auto res = matrix_vector_multiply(matrix, float4(pt[0], pt[1], pt[2], 1.f));
    return vec3<float>(res[0], res[1], res[2]);
    }

  JTKQBVHDEF vec3<float> transform_vector(const float4x4& matrix, const vec3<float>& vec)
    {
    auto res = matrix_vector_multiply(matrix, float4(vec[0], vec[1], vec[2], 0.f));
    return vec3<float>(res[0], res[1], res[2]);
    }

  JTKQBVHDEF vec3<float> transform(const float4x4& matrix, const vec3<float>& pt, bool is_vector)
    {
    if (is_vector)
      return transform_vector(matrix, pt);
    else
      return transform(matrix, pt);
    }

  JTKQBVHDEF float4 transform(const float4x4& matrix, const float4& pt)
    {
    auto res = matrix_vector_multiply(matrix, pt);
    if (res[3] != 1.f && res[3])
      {
      res[0] /= res[3];
      res[1] /= res[3];
      res[2] /= res[3];
      res[3] = 1.f;
      }
    return res;
    }

  JTKQBVHDEF float4x4 make_transformation(const vec3<float>& i_origin, const vec3<float>& i_x_axis, const vec3<float>& i_y_axis, const vec3<float>& i_z_axis)
    {
    float4x4 matrix;
    matrix[0] = i_x_axis[0];
    matrix[1] = i_x_axis[1];
    matrix[2] = i_x_axis[2];
    matrix[3] = 0;
    matrix[4] = i_y_axis[0];
    matrix[5] = i_y_axis[1];
    matrix[6] = i_y_axis[2];
    matrix[7] = 0;
    matrix[8] = i_z_axis[0];
    matrix[9] = i_z_axis[1];
    matrix[10] = i_z_axis[2];
    matrix[11] = 0;
    matrix[12] = i_origin[0];
    matrix[13] = i_origin[1];
    matrix[14] = i_origin[2];
    matrix[15] = 1;
    return matrix;
    }

  JTKQBVHDEF float4x4 make_rotation(const vec3<float>& i_position, const vec3<float>& i_direction, float i_angle_radians)
    {
    auto matrix = make_identity();
    auto direction = normalize(i_direction);

    auto cos_alpha = std::cos(i_angle_radians);
    auto sin_alpha = std::sin(i_angle_radians);

    matrix[0] = (direction[0] * direction[0]) * (1 - cos_alpha) + cos_alpha;
    matrix[4] = (direction[0] * direction[1]) * (1 - cos_alpha) - direction[2] * sin_alpha;
    matrix[8] = (direction[0] * direction[2]) * (1 - cos_alpha) + direction[1] * sin_alpha;

    matrix[1] = (direction[0] * direction[1]) * (1 - cos_alpha) + direction[2] * sin_alpha;
    matrix[5] = (direction[1] * direction[1]) * (1 - cos_alpha) + cos_alpha;
    matrix[9] = (direction[1] * direction[2]) * (1 - cos_alpha) - direction[0] * sin_alpha;

    matrix[2] = (direction[0] * direction[2]) * (1 - cos_alpha) - direction[1] * sin_alpha;
    matrix[6] = (direction[1] * direction[2]) * (1 - cos_alpha) + direction[0] * sin_alpha;
    matrix[10] = (direction[2] * direction[2]) * (1 - cos_alpha) + cos_alpha;

    auto rotated_position = transform_vector(matrix, i_position);

    matrix[12] = i_position[0] - rotated_position[0];
    matrix[13] = i_position[1] - rotated_position[1];
    matrix[14] = i_position[2] - rotated_position[2];

    return matrix;
    }

  JTKQBVHDEF float4x4 make_scale3d(float scale_x, float scale_y, float scale_z)
    {
    return float4x4(float4(scale_x, 0.f, 0.f, 0.f), float4(0.f, scale_y, 0.f, 0.f), float4(0.f, 0.f, scale_z, 0.f), float4(0.f, 0.f, 0.f, 1.f));
    }

  JTKQBVHDEF float4x4 make_translation(const vec3<float>& i_translation)
    {
    return make_translation(i_translation[0], i_translation[1], i_translation[2]);
    }

  JTKQBVHDEF vec3<float> get_translation(const float4x4& matrix)
    {
    return vec3<float>(matrix[12], matrix[13], matrix[14]);
    }

  JTKQBVHDEF void set_x_axis(float4x4& matrix, const vec3<float>& x)
    {
    matrix[0] = x[0];
    matrix[1] = x[1];
    matrix[2] = x[2];
    }

  JTKQBVHDEF void set_y_axis(float4x4& matrix, const vec3<float>& y)
    {
    matrix[4] = y[0];
    matrix[5] = y[1];
    matrix[6] = y[2];
    }

  JTKQBVHDEF void set_z_axis(float4x4& matrix, const vec3<float>& z)
    {
    matrix[8] = z[0];
    matrix[9] = z[1];
    matrix[10] = z[2];
    }

  JTKQBVHDEF void set_translation(float4x4& matrix, const vec3<float>& t)
    {
    matrix[12] = t[0];
    matrix[13] = t[1];
    matrix[14] = t[2];
    }

  JTKQBVHDEF vec3<float> get_x_axis(const float4x4& matrix)
    {
    return vec3<float>(matrix[0], matrix[1], matrix[2]);
    }

  JTKQBVHDEF vec3<float> get_y_axis(const float4x4& matrix)
    {
    return vec3<float>(matrix[4], matrix[5], matrix[6]);
    }

  JTKQBVHDEF vec3<float> get_z_axis(const float4x4& matrix)
    {
    return vec3<float>(matrix[8], matrix[9], matrix[10]);
    }

  JTKQBVHDEF float determinant(const float4x4& m)
    {
    auto inv0 = m[5] * m[10] * m[15] -
      m[5] * m[11] * m[14] -
      m[9] * m[6] * m[15] +
      m[9] * m[7] * m[14] +
      m[13] * m[6] * m[11] -
      m[13] * m[7] * m[10];

    auto inv4 = -m[4] * m[10] * m[15] +
      m[4] * m[11] * m[14] +
      m[8] * m[6] * m[15] -
      m[8] * m[7] * m[14] -
      m[12] * m[6] * m[11] +
      m[12] * m[7] * m[10];

    auto inv8 = m[4] * m[9] * m[15] -
      m[4] * m[11] * m[13] -
      m[8] * m[5] * m[15] +
      m[8] * m[7] * m[13] +
      m[12] * m[5] * m[11] -
      m[12] * m[7] * m[9];

    auto inv12 = -m[4] * m[9] * m[14] +
      m[4] * m[10] * m[13] +
      m[8] * m[5] * m[14] -
      m[8] * m[6] * m[13] -
      m[12] * m[5] * m[10] +
      m[12] * m[6] * m[9];

    return m[0] * inv0 + m[1] * inv4 + m[2] * inv8 + m[3] * inv12;
    }

  /////////////////////////////////////////////////////////////////////////
  // qbvh
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF qbvh_voxel* build_triangle_qbvh_voxels(qbvh_voxel& total_bb, qbvh_voxel& centroid_bb, const vec3<float>* vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
    {
    qbvh_voxel* lst = new qbvh_voxel[nr_of_triangles];
    total_bb.bbox_min = std::numeric_limits<float>::max();
    total_bb.bbox_max = -std::numeric_limits<float>::max();
    centroid_bb.bbox_min = total_bb.bbox_min;
    centroid_bb.bbox_max = total_bb.bbox_max;
    const unsigned int number_of_threads = hardware_concurrency();

    aligned_vector<qbvh_voxel> total_bbs(number_of_threads, total_bb);
    aligned_vector<qbvh_voxel> total_centroids(number_of_threads, centroid_bb);
    parallel_for((unsigned int)0, number_of_threads, [&](unsigned int i)
      {
      uint32_t s = (uint32_t)((uint64_t)i * (uint64_t)(nr_of_triangles) / (uint64_t)number_of_threads);
      const uint32_t e = (uint32_t)((uint64_t)(i + 1) * (uint64_t)(nr_of_triangles) / (uint64_t)number_of_threads);
      for (uint32_t t = s; t < e; ++t)
        {
        const float4 v0(vertices[triangles[t][0]][0], vertices[triangles[t][0]][1], vertices[triangles[t][0]][2], 1);
        const float4 v1(vertices[triangles[t][1]][0], vertices[triangles[t][1]][1], vertices[triangles[t][1]][2], 1);
        const float4 v2(vertices[triangles[t][2]][0], vertices[triangles[t][2]][1], vertices[triangles[t][2]][2], 1);
        lst[t].bbox_min = min(min(v0, v1), v2);
        lst[t].bbox_max = max(max(v0, v1), v2);
        lst[t].centroid = (lst[t].bbox_min + lst[t].bbox_max)*0.5;
        total_bbs[i].bbox_min = min(total_bbs[i].bbox_min, lst[t].bbox_min);
        total_bbs[i].bbox_max = max(total_bbs[i].bbox_max, lst[t].bbox_max);
        total_centroids[i].bbox_min = min(total_centroids[i].bbox_min, lst[t].centroid);
        total_centroids[i].bbox_max = max(total_centroids[i].bbox_max, lst[t].centroid);
        }
      });

    for (uint32_t i = 0; i < number_of_threads; ++i)
      {
      total_bb.bbox_min = min(total_bb.bbox_min, total_bbs[i].bbox_min);
      total_bb.bbox_max = max(total_bb.bbox_max, total_bbs[i].bbox_max);
      centroid_bb.bbox_min = min(centroid_bb.bbox_min, total_centroids[i].bbox_min);
      centroid_bb.bbox_max = max(centroid_bb.bbox_max, total_centroids[i].bbox_max);
      }
    total_bb.centroid = (total_bb.bbox_min + total_bb.bbox_max)*0.5;
    centroid_bb.centroid = (centroid_bb.bbox_min + centroid_bb.bbox_max)*0.5;
    return lst;
    }

  JTKQBVHDEF qbvh_voxel* build_sphere_qbvh_voxels(qbvh_voxel& total_bb, qbvh_voxel& centroid_bb, const vec3<float>* origins, const float* radii, uint32_t nr_of_spheres)
    {
    qbvh_voxel* lst = new qbvh_voxel[nr_of_spheres];
    total_bb.bbox_min = std::numeric_limits<float>::max();
    total_bb.bbox_max = -std::numeric_limits<float>::max();
    centroid_bb.bbox_min = total_bb.bbox_min;
    centroid_bb.bbox_max = total_bb.bbox_max;
    const unsigned int number_of_threads = hardware_concurrency();

    aligned_vector<qbvh_voxel> total_bbs(number_of_threads, total_bb);
    aligned_vector<qbvh_voxel> total_centroids(number_of_threads, centroid_bb);
    parallel_for((unsigned int)0, number_of_threads, [&](unsigned int i)
      {
      uint32_t s = (uint32_t)((uint64_t)i * (uint64_t)(nr_of_spheres) / (uint64_t)number_of_threads);
      const uint32_t e = (uint32_t)((uint64_t)(i + 1) * (uint64_t)(nr_of_spheres) / (uint64_t)number_of_threads);
      for (uint32_t t = s; t < e; ++t)
        {
        const float4 v(origins[t][0], origins[t][1], origins[t][2], 1);
        lst[t].bbox_min = v - float4(radii[t], radii[t], radii[t], 0.f);
        lst[t].bbox_max = v + float4(radii[t], radii[t], radii[t], 0.f);
        lst[t].centroid = v;
        total_bbs[i].bbox_min = min(total_bbs[i].bbox_min, lst[t].bbox_min);
        total_bbs[i].bbox_max = max(total_bbs[i].bbox_max, lst[t].bbox_max);
        total_centroids[i].bbox_min = min(total_centroids[i].bbox_min, lst[t].centroid);
        total_centroids[i].bbox_max = max(total_centroids[i].bbox_max, lst[t].centroid);
        }
      });

    for (uint32_t i = 0; i < number_of_threads; ++i)
      {
      total_bb.bbox_min = min(total_bb.bbox_min, total_bbs[i].bbox_min);
      total_bb.bbox_max = max(total_bb.bbox_max, total_bbs[i].bbox_max);
      centroid_bb.bbox_min = min(centroid_bb.bbox_min, total_centroids[i].bbox_min);
      centroid_bb.bbox_max = max(centroid_bb.bbox_max, total_centroids[i].bbox_max);
      }
    total_bb.centroid = (total_bb.bbox_min + total_bb.bbox_max)*0.5;
    centroid_bb.centroid = (centroid_bb.bbox_min + centroid_bb.bbox_max)*0.5;
    return lst;
    }

  JTKQBVHDEF void unite_four_aabbs(float4* out, float4* in, int k)
    {
    out[0][k] = min_horizontal(in[0]);
    out[1][k] = min_horizontal(in[1]);
    out[2][k] = min_horizontal(in[2]);
    out[3][k] = max_horizontal(in[3]);
    out[4][k] = max_horizontal(in[4]);
    out[5][k] = max_horizontal(in[5]);
    }


  JTKQBVHDEF void get_bbox(float4* bbox, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last, int k)
    {
    assert(first != last);
    bbox[0][k] = voxels[*first].bbox_min[0];
    bbox[1][k] = voxels[*first].bbox_min[1];
    bbox[2][k] = voxels[*first].bbox_min[2];
    bbox[3][k] = voxels[*first].bbox_max[0];
    bbox[4][k] = voxels[*first].bbox_max[1];
    bbox[5][k] = voxels[*first].bbox_max[2];
    for (auto it = first + 1; it != last; ++it)
      {
      bbox[0][k] = std::min<float>(bbox[0][k], voxels[*it].bbox_min[0]);
      bbox[1][k] = std::min<float>(bbox[1][k], voxels[*it].bbox_min[1]);
      bbox[2][k] = std::min<float>(bbox[2][k], voxels[*it].bbox_min[2]);
      bbox[3][k] = std::max<float>(bbox[3][k], voxels[*it].bbox_max[0]);
      bbox[4][k] = std::max<float>(bbox[4][k], voxels[*it].bbox_max[1]);
      bbox[5][k] = std::max<float>(bbox[5][k], voxels[*it].bbox_max[2]);
      }
    }


  JTKQBVHDEF float calculate_half_surface_area(const float4& left, const float4& right)
    {
    const float4 diff = right - left;
    const float4 diff2(diff[1], diff[2], diff[0], 0);
    float4 tmp = diff * diff2;
    return tmp[0] + tmp[1] + tmp[2];
    }


  JTKQBVHDEF uint8_t find_largest_dimension(const qbvh_voxel& bb)
    {
    const float4 diff = bb.bbox_max - bb.bbox_min;
    const uint8_t dim = diff[1] > diff[0] ? (diff[2] > diff[1] ? 2 : 1) : (diff[2] > diff[0] ? 2 : 0);
    return dim;
    }


  JTKQBVHDEF std::vector<uint32_t>::iterator partition(qbvh_voxel& bbox_left, qbvh_voxel& bbox_right, qbvh_voxel& centroid_left, qbvh_voxel& centroid_right, const uint32_t dim, const float split_pos, const qbvh_voxel* voxels, std::vector<uint32_t>::iterator first, std::vector<uint32_t>::iterator last)
    {
    auto last_save = last;
    while (first != last)
      {
      while (voxels[*first].centroid[dim] < split_pos)
        {
        bbox_left.bbox_min = min(bbox_left.bbox_min, voxels[*first].bbox_min);
        bbox_left.bbox_max = max(bbox_left.bbox_max, voxels[*first].bbox_max);
        centroid_left.bbox_min = min(centroid_left.bbox_min, voxels[*first].centroid);
        centroid_left.bbox_max = max(centroid_left.bbox_max, voxels[*first].centroid);
        ++first;
        if (first == last)
          return first;
        }
      do
        {
        if (last != last_save)
          {
          bbox_right.bbox_min = min(bbox_right.bbox_min, voxels[*last].bbox_min);
          bbox_right.bbox_max = max(bbox_right.bbox_max, voxels[*last].bbox_max);
          centroid_right.bbox_min = min(centroid_right.bbox_min, voxels[*last].centroid);
          centroid_right.bbox_max = max(centroid_right.bbox_max, voxels[*last].centroid);
          }
        --last;
        if (first == last)
          {
          bbox_right.bbox_min = min(bbox_right.bbox_min, voxels[*last].bbox_min);
          bbox_right.bbox_max = max(bbox_right.bbox_max, voxels[*last].bbox_max);
          centroid_right.bbox_min = min(centroid_right.bbox_min, voxels[*last].centroid);
          centroid_right.bbox_max = max(centroid_right.bbox_max, voxels[*last].centroid);
          return first;
          }
        } while (voxels[*last].centroid[dim] >= split_pos);
        std::swap(*first, *last);
        bbox_left.bbox_min = min(bbox_left.bbox_min, voxels[*first].bbox_min);
        bbox_left.bbox_max = max(bbox_left.bbox_max, voxels[*first].bbox_max);
        centroid_left.bbox_min = min(centroid_left.bbox_min, voxels[*first].centroid);
        centroid_left.bbox_max = max(centroid_left.bbox_max, voxels[*first].centroid);
        bbox_right.bbox_min = min(bbox_right.bbox_min, voxels[*last].bbox_min);
        bbox_right.bbox_max = max(bbox_right.bbox_max, voxels[*last].bbox_max);
        centroid_right.bbox_min = min(centroid_right.bbox_min, voxels[*last].centroid);
        centroid_right.bbox_max = max(centroid_right.bbox_max, voxels[*last].centroid);
        ++first;
      }
    return first;
    }

  /////////////////////////////////////////////////////////////////////////
  // qbvh_two_level
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF void unite_four_aabbs(qbvh_voxel& out, const float4* in)
    {
    out.bbox_min[0] = min_horizontal(in[0]);
    out.bbox_min[1] = min_horizontal(in[1]);
    out.bbox_min[2] = min_horizontal(in[2]);
    out.bbox_min[3] = 1.f;
    out.bbox_max[0] = max_horizontal(in[3]);
    out.bbox_max[1] = max_horizontal(in[4]);
    out.bbox_max[2] = max_horizontal(in[5]);
    out.bbox_max[3] = 1.f;
    }

  /////////////////////////////////////////////////////////////////////////
  // qbvh_two_level_with_transformations
  /////////////////////////////////////////////////////////////////////////

  JTKQBVHDEF void transform(qbvh_voxel& voxel, const float4x4& matrix)
    {
    auto xa = matrix.col[0] * voxel.bbox_min[0];
    auto xb = matrix.col[0] * voxel.bbox_max[0];

    auto ya = matrix.col[1] * voxel.bbox_min[1];
    auto yb = matrix.col[1] * voxel.bbox_max[1];

    auto za = matrix.col[2] * voxel.bbox_min[2];
    auto zb = matrix.col[2] * voxel.bbox_max[2];

    voxel.bbox_min = min(xa, xb) + min(ya, yb) + min(za, zb) + matrix.col[3];
    voxel.bbox_max = max(xa, xb) + max(ya, yb) + max(za, zb) + matrix.col[3];
    }


  } // namespace jtk

#endif //JTK_QBVH_IMPLEMENTATION
