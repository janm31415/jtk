#pragma once

#include "concurrency.h"

#include <immintrin.h>

#include <vector>

//#define _AVX2

// software rendering of point clouds

/*
TO CHECK: std::vector is not necessarily aligned. Might need to replace std::vector by similar type with aligned memory
*/

namespace jtk
  {

  typedef union
    {
    struct
      {
      unsigned char red;
      unsigned char green;
      unsigned char blue;
      unsigned char alpha;
      } rgba;
    uint32_t color;
    } rgb_value;

  inline uint32_t make_color(unsigned char alpha, unsigned char red, unsigned char green, unsigned char blue)
    {
    return (uint32_t(alpha) << 24) | (uint32_t(blue) << 16) | (uint32_t(green) << 8) | uint32_t(red);
    }

  inline uint32_t make_color(unsigned char red, unsigned char green, unsigned char blue)
    {
    return 0xff000000 | (uint32_t(blue) << 16) | (uint32_t(green) << 8) | uint32_t(red);
    }

  inline uint32_t get_red(uint32_t color)
    {
    return (color & 0x000000ff);
    }

  inline uint32_t get_green(uint32_t color)
    {
    return (color & 0x0000ff00) >> 8;
    }

  inline uint32_t get_blue(uint32_t color)
    {
    return (color & 0x00ff0000) >> 16;
    }

  inline uint32_t get_alpha(uint32_t color)
    {
    return (color & 0xff000000) >> 24;
    }

  struct frame_buffer
    {
    int w;
    int h;
    uint32_t* pixels;
    float* zbuffer;
    };

  struct object_buffer
    {
    object_buffer() : number_of_vertices(0) {}
    const float* vertices;
    const float* normals;
    const uint32_t* colors;
    uint32_t number_of_vertices;
    };

#ifdef _WIN32
  __declspec(align(32))
#endif
  struct render_data
    {
#ifdef _AVX2
    __m256 projection_times_camera_times_object_avx[16];
#else
    __m128 projection_times_camera_times_object_sse[16];
#endif
    frame_buffer fb;
    object_buffer obj;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<uint32_t> vertex_clip_info;
    float camera_position[16];
    float object_system[16];
    float projection_matrix[16];
    float projection_times_camera_times_object[16];
    uint32_t ob_id;
    float light_direction[3];
    uint32_t light_color;
    int point_size;
    uint16_t pad1;
    render_data()
      {
      point_size = 0;
      ob_id = 0;
      light_direction[0] = 0.f;
      light_direction[1] = 0.f;
      light_direction[2] = 1.f;
      light_color = 0xffffffff;
      }
    }
#ifndef _WIN32 // linux alignment in gcc
  __attribute__((aligned(32)))
#endif
    ;

  inline float get_depth(const render_data& rd, uint32_t x, uint32_t y)
    {
    uint32_t index = x + y * rd.fb.w;
    float z = 1.f / rd.fb.zbuffer[index];
    return z;
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

  inline void invert_projection_matrix(float* out, const float* in)
    {
    out[0] = 1.f / in[0];
    out[1] = 0.f;
    out[2] = 0.f;
    out[3] = 0.f;
    out[4] = 0.f;
    out[5] = 1.f / in[5];
    out[6] = 0.f;
    out[7] = 0.f;
    out[8] = 0.f;
    out[9] = 0.f;
    out[10] = 0.f;
    out[11] = 1.f / in[14];
    out[12] = in[8] / in[0];
    out[13] = in[9] / in[5];
    out[14] = -1.f;
    out[15] = in[10] / in[14];
    }

  inline void frustum(float* out, float left, float right, float bottom, float top, float _near, float _far)
    {
    out[0 + 4 * 0] = 2.f*_near / (right - left);
    out[0 + 4 * 1] = 0.f;
    out[0 + 4 * 2] = (right + left) / (right - left);
    out[0 + 4 * 3] = 0.f;
    out[1 + 4 * 0] = 0.f;
    out[1 + 4 * 1] = -2.f*_near / (top - bottom);
    out[1 + 4 * 2] = -(top + bottom) / (top - bottom);
    out[1 + 4 * 3] = 0.f;
    out[2 + 4 * 0] = 0.f;
    out[2 + 4 * 1] = 0.f;
    out[2 + 4 * 2] = -(_far + _near) / (_far - _near);
    out[2 + 4 * 3] = -(2.f * _far * _near) / (_far - _near);
    out[3 + 4 * 0] = 0.f;
    out[3 + 4 * 1] = 0.f;
    out[3 + 4 * 2] = -1.f;
    out[3 + 4 * 3] = 0.f;
    }

  inline void get_identity(float* m)
    {
    m[0] = 1.f;
    m[1] = 0.f;
    m[2] = 0.f;
    m[3] = 0.f;
    m[4] = 0.f;
    m[5] = 1.f;
    m[6] = 0.f;
    m[7] = 0.f;
    m[8] = 0.f;
    m[9] = 0.f;
    m[10] = 1.f;
    m[11] = 0.f;
    m[12] = 0.f;
    m[13] = 0.f;
    m[14] = 0.f;
    m[15] = 1.f;
    }

  inline void make_translation(float* m, float x, float y, float z)
    {
    get_identity(m);
    m[12] = x;
    m[13] = y;
    m[14] = z;
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

  inline void mat_vec_multiply(float* out, const float* m, const float* v)
    {
    out[0] = m[0] * v[0] + m[4] * v[1] + m[8] * v[2] + m[12] * v[3];
    out[1] = m[1] * v[0] + m[5] * v[1] + m[9] * v[2] + m[13] * v[3];
    out[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3];
    out[3] = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
    }

  inline void clear_frame_buffer(render_data& rd, uint32_t rgba, float depth)
    {
    std::fill(rd.fb.pixels, rd.fb.pixels + rd.fb.w*rd.fb.h, rgba);
    std::fill(rd.fb.zbuffer, rd.fb.zbuffer + rd.fb.w*rd.fb.h, 1.f / depth);
    }

  inline float min2(float a, float b)
    {
    return a < b ? a : b;
    }

  inline void bind(render_data& rd, uint32_t object_buffer_id)
    {
    rd.ob_id = object_buffer_id;
    }

  inline void bind(render_data& rd, const frame_buffer& fb)
    {
    rd.fb = fb;
    }

  inline void bind(render_data& rd, const object_buffer& ob)
    {
    rd.obj = ob;
    rd.vertices_x.reserve(ob.number_of_vertices);
    rd.vertices_y.reserve(ob.number_of_vertices);
    rd.vertices_z.reserve(ob.number_of_vertices);
    rd.vertex_clip_info.reserve(ob.number_of_vertices);
    }

  inline void bind(render_data& rd, const float* camera_position, const float* object_system, const float* projection_matrix)
    {
    std::memcpy(rd.camera_position, camera_position, 16 * sizeof(float));
    std::memcpy(rd.object_system, object_system, 16 * sizeof(float));
    std::memcpy(rd.projection_matrix, projection_matrix, 16 * sizeof(float));
    float temp[16];
    matrix_multiply(temp, camera_position, object_system);
    matrix_multiply(rd.projection_times_camera_times_object, projection_matrix, temp);
    for (size_t i = 0; i < 16; ++i)
#ifdef AVX2
      rd.projection_times_camera_times_object_avx[i] = _mm256_set1_ps(rd.projection_times_camera_times_object[i]);
#else
      rd.projection_times_camera_times_object_sse[i] = _mm_set1_ps(rd.projection_times_camera_times_object[i]);
#endif

    }

  inline void world_to_screen(float& x, float& y, float& z, const render_data& rd)
    {
    float v[4], v2[4];
    v[0] = x;
    v[1] = y;
    v[2] = z;
    v[3] = 1.f;
    mat_vec_multiply(v2, rd.camera_position, v);
    mat_vec_multiply(v, rd.projection_matrix, v2);
    v[0] /= v[3];
    v[1] /= v[3];
    v[2] /= v[3];
    x = (v[0] + 1.f)*rd.fb.w * 0.5f;
    y = (v[1] + 1.f)*rd.fb.h * 0.5f;
    z = v[3];
    }

  inline void _draw(render_data& rd)
    {
    rd.vertices_x.clear();
    rd.vertices_y.clear();
    rd.vertices_z.clear();
    rd.vertex_clip_info.clear();
    const float* p_vert = rd.obj.vertices;
#ifdef _AVX2
    uint32_t sz = rd.obj.number_of_vertices - (rd.obj.number_of_vertices & 7);
#else
    uint32_t sz = rd.obj.number_of_vertices - (rd.obj.number_of_vertices & 3);
#endif
    rd.vertices_x.resize(rd.obj.number_of_vertices);
    rd.vertices_y.resize(rd.obj.number_of_vertices);
    rd.vertices_z.resize(rd.obj.number_of_vertices);
    rd.vertex_clip_info.resize(rd.obj.number_of_vertices);
    //for (uint32_t i = 0; i < sz; i += 8)

    const unsigned int number_of_threads = hardware_concurrency() * 4;

    parallel_for((unsigned int)0, number_of_threads, [&](unsigned int i)
      {
#ifdef AVX2
      const uint32_t s = (uint32_t)((uint64_t)i * (uint64_t)(sz / 8) / (uint64_t)number_of_threads);
      const uint32_t e = (uint32_t)((uint64_t)(i + 1) * (uint64_t)(sz / 8) / (uint64_t)number_of_threads);
      for (uint32_t t = s; t < e; ++t)
        {
        //https://software.intel.com/content/www/us/en/develop/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx.html
        //Array of Structs to Struct Of Arrays

        const float* p = p_vert + t * 8 * 3;
        const __m128 *m = (const __m128*) p;
        __m256 m03;
        __m256 m14;
        __m256 m25;
        m03 = _mm256_castps128_ps256(m[0]); // load lower halves
        m14 = _mm256_castps128_ps256(m[1]);
        m25 = _mm256_castps128_ps256(m[2]);
        m03 = _mm256_insertf128_ps(m03, m[3], 1);  // load upper halves
        m14 = _mm256_insertf128_ps(m14, m[4], 1);
        m25 = _mm256_insertf128_ps(m25, m[5], 1);

        __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE(2, 1, 3, 2)); // upper x's and y's 
        __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE(1, 0, 2, 1)); // lower y's and z's
        __m256 X = _mm256_shuffle_ps(m03, xy, _MM_SHUFFLE(2, 0, 3, 0));
        __m256 Y = _mm256_shuffle_ps(yz, xy, _MM_SHUFFLE(3, 1, 2, 0));
        __m256 Z = _mm256_shuffle_ps(yz, m25, _MM_SHUFFLE(3, 0, 3, 1));
        __m256 W = _mm256_set1_ps(1.f);

        //p_vert += 3 * 8;


        __m256 VX = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rd.projection_times_camera_times_object_avx[0], X), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[4], Y)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[8], Z)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[12], W));
        __m256 VY = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rd.projection_times_camera_times_object_avx[1], X), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[5], Y)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[9], Z)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[13], W));
        __m256 VZ = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rd.projection_times_camera_times_object_avx[2], X), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[6], Y)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[10], Z)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[14], W));
        __m256 VW = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rd.projection_times_camera_times_object_avx[3], X), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[7], Y)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[11], Z)), _mm256_mul_ps(rd.projection_times_camera_times_object_avx[15], W));
        __m256i clip = _mm256_set1_epi32(0);
        VX = _mm256_div_ps(VX, VW);
        VY = _mm256_div_ps(VY, VW);
        VZ = _mm256_div_ps(VZ, VW);
        __m256i mxlt = _mm256_castps_si256(_mm256_cmp_ps(VX, _mm256_set1_ps(-1.f), 1));
        __m256i mxgt = _mm256_castps_si256(_mm256_cmp_ps(VX, _mm256_set1_ps(1.f), 14));
        __m256i mylt = _mm256_castps_si256(_mm256_cmp_ps(VY, _mm256_set1_ps(-1.f), 1));
        __m256i mygt = _mm256_castps_si256(_mm256_cmp_ps(VY, _mm256_set1_ps(1.f), 14));
        __m256i mzlt = _mm256_castps_si256(_mm256_cmp_ps(VZ, _mm256_set1_ps(-1.f), 1));
        __m256i mzgt = _mm256_castps_si256(_mm256_cmp_ps(VZ, _mm256_set1_ps(1.f), 14));

        clip = _mm256_or_si256(clip, _mm256_and_si256(mxlt, _mm256_set1_epi32(1)));
        clip = _mm256_or_si256(clip, _mm256_and_si256(mxgt, _mm256_set1_epi32(2)));
        clip = _mm256_or_si256(clip, _mm256_and_si256(mylt, _mm256_set1_epi32(4)));
        clip = _mm256_or_si256(clip, _mm256_and_si256(mygt, _mm256_set1_epi32(8)));
        clip = _mm256_or_si256(clip, _mm256_and_si256(mzlt, _mm256_set1_epi32(16)));
        clip = _mm256_or_si256(clip, _mm256_and_si256(mzgt, _mm256_set1_epi32(32)));

        VX = _mm256_mul_ps(_mm256_add_ps(VX, _mm256_set1_ps(1.f)), _mm256_set1_ps(rd.fb.w * 0.5f));
        VY = _mm256_mul_ps(_mm256_add_ps(VY, _mm256_set1_ps(1.f)), _mm256_set1_ps(rd.fb.h * 0.5f));

        const float* vx = rd.vertices_x.data() + t * 8;
        __m256* m2 = (__m256*)vx;
        m2[0] = VX;
        const float* vy = rd.vertices_y.data() + t * 8;
        m2 = (__m256*)vy;
        m2[0] = VY;
        const float* vz = rd.vertices_z.data() + t * 8;
        m2 = (__m256*)vz;
        m2[0] = VW;

        uint32_t* c = rd.vertex_clip_info.data() + t * 8;
        __m256i* mm = (__m256i*)c;
        mm[0] = clip;
        }
#else
      const uint32_t s = (uint32_t)((uint64_t)i * (uint64_t)(sz / 4) / (uint64_t)number_of_threads);
      const uint32_t e = (uint32_t)((uint64_t)(i + 1) * (uint64_t)(sz / 4) / (uint64_t)number_of_threads);
      for (uint32_t t = s; t < e; ++t)
        {
        //https://software.intel.com/content/www/us/en/develop/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx.html
        //Array of Structs to Struct Of Arrays


        /*
        float *p;  // address of first vector
 __m128 x0y0z0x1 = _mm_load_ps(p+0);
 __m128 y1z1x2y2 = _mm_load_ps(p+4);
 __m128 z2x3y3z3 = _mm_load_ps(p+8);
 __m128 x2y2x3y3 = _mm_shuffle_ps(y1z1x2y2,z2x3y3z3,_MM_SHUFFLE( 2,1,3,2));
 __m128 y0z0y1z1 = _mm_shuffle_ps(x0y0z0x1,y1z1x2y2,_MM_SHUFFLE( 1,0,2,1));
 __m128 x        = _mm_shuffle_ps(x0y0z0x1,x2y2x3y3,_MM_SHUFFLE( 2,0,3,0)); // x0x1x2x3
 __m128 y        = _mm_shuffle_ps(y0z0y1z1,x2y2x3y3,_MM_SHUFFLE( 3,1,2,0)); // y0y1y2y3
 __m128 z        = _mm_shuffle_ps(y0z0y1z1,z2x3y3z3,_MM_SHUFFLE( 3,0,3,1)); // z0z1z2z3
        */

        const float* p = p_vert + t * 4 * 3;
        __m128 x0y0z0x1 = _mm_load_ps(p + 0);
        __m128 y1z1x2y2 = _mm_load_ps(p + 4);
        __m128 z2x3y3z3 = _mm_load_ps(p + 8);
        __m128 x2y2x3y3 = _mm_shuffle_ps(y1z1x2y2, z2x3y3z3, _MM_SHUFFLE(2, 1, 3, 2));
        __m128 y0z0y1z1 = _mm_shuffle_ps(x0y0z0x1, y1z1x2y2, _MM_SHUFFLE(1, 0, 2, 1));
        __m128 X = _mm_shuffle_ps(x0y0z0x1, x2y2x3y3, _MM_SHUFFLE(2, 0, 3, 0)); // x0x1x2x3
        __m128 Y = _mm_shuffle_ps(y0z0y1z1, x2y2x3y3, _MM_SHUFFLE(3, 1, 2, 0)); // y0y1y2y3
        __m128 Z = _mm_shuffle_ps(y0z0y1z1, z2x3y3z3, _MM_SHUFFLE(3, 0, 3, 1)); // z0z1z2z3
        __m128 W = _mm_set1_ps(1.f);


        __m128 VX = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(rd.projection_times_camera_times_object_sse[0], X), _mm_mul_ps(rd.projection_times_camera_times_object_sse[4], Y)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[8], Z)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[12], W));
        __m128 VY = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(rd.projection_times_camera_times_object_sse[1], X), _mm_mul_ps(rd.projection_times_camera_times_object_sse[5], Y)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[9], Z)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[13], W));
        __m128 VZ = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(rd.projection_times_camera_times_object_sse[2], X), _mm_mul_ps(rd.projection_times_camera_times_object_sse[6], Y)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[10], Z)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[14], W));
        __m128 VW = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(rd.projection_times_camera_times_object_sse[3], X), _mm_mul_ps(rd.projection_times_camera_times_object_sse[7], Y)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[11], Z)), _mm_mul_ps(rd.projection_times_camera_times_object_sse[15], W));
        __m128i clip = _mm_set1_epi32(0);
        VX = _mm_div_ps(VX, VW);
        VY = _mm_div_ps(VY, VW);
        VZ = _mm_div_ps(VZ, VW);
        __m128i mxlt = _mm_castps_si128(_mm_cmplt_ps(VX, _mm_set1_ps(-1.f)));
        __m128i mxgt = _mm_castps_si128(_mm_cmpgt_ps(VX, _mm_set1_ps(1.f)));
        __m128i mylt = _mm_castps_si128(_mm_cmplt_ps(VY, _mm_set1_ps(-1.f)));
        __m128i mygt = _mm_castps_si128(_mm_cmpgt_ps(VY, _mm_set1_ps(1.f)));
        __m128i mzlt = _mm_castps_si128(_mm_cmplt_ps(VZ, _mm_set1_ps(-1.f)));
        __m128i mzgt = _mm_castps_si128(_mm_cmpgt_ps(VZ, _mm_set1_ps(1.f)));

        clip = _mm_or_si128(clip, _mm_and_si128(mxlt, _mm_set1_epi32(1)));
        clip = _mm_or_si128(clip, _mm_and_si128(mxgt, _mm_set1_epi32(2)));
        clip = _mm_or_si128(clip, _mm_and_si128(mylt, _mm_set1_epi32(4)));
        clip = _mm_or_si128(clip, _mm_and_si128(mygt, _mm_set1_epi32(8)));
        clip = _mm_or_si128(clip, _mm_and_si128(mzlt, _mm_set1_epi32(16)));
        clip = _mm_or_si128(clip, _mm_and_si128(mzgt, _mm_set1_epi32(32)));

        VX = _mm_mul_ps(_mm_add_ps(VX, _mm_set1_ps(1.f)), _mm_set1_ps(rd.fb.w * 0.5f));
        VY = _mm_mul_ps(_mm_add_ps(VY, _mm_set1_ps(1.f)), _mm_set1_ps(rd.fb.h * 0.5f));

        const float* vx = rd.vertices_x.data() + t * 4;
        __m128* m2 = (__m128*)vx;
        m2[0] = VX;
        const float* vy = rd.vertices_y.data() + t * 4;
        m2 = (__m128*)vy;
        m2[0] = VY;
        const float* vz = rd.vertices_z.data() + t * 4;
        m2 = (__m128*)vz;
        m2[0] = VW;

        uint32_t* c = rd.vertex_clip_info.data() + t * 4;
        __m128i* mm = (__m128i*)c;
        mm[0] = clip;
        }
#endif
      });

    for (uint32_t i = sz; i < rd.obj.number_of_vertices; ++i)
      {

      float X = *p_vert++;
      float Y = *p_vert++;
      float Z = *p_vert++;

      float VX = rd.projection_times_camera_times_object[0] * X + rd.projection_times_camera_times_object[4] * Y + rd.projection_times_camera_times_object[8] * Z + rd.projection_times_camera_times_object[12];
      float VY = rd.projection_times_camera_times_object[1] * X + rd.projection_times_camera_times_object[5] * Y + rd.projection_times_camera_times_object[9] * Z + rd.projection_times_camera_times_object[13];
      float VZ = rd.projection_times_camera_times_object[2] * X + rd.projection_times_camera_times_object[6] * Y + rd.projection_times_camera_times_object[10] * Z + rd.projection_times_camera_times_object[14];
      float VW = rd.projection_times_camera_times_object[3] * X + rd.projection_times_camera_times_object[7] * Y + rd.projection_times_camera_times_object[11] * Z + rd.projection_times_camera_times_object[15];
      uint32_t clip_info = 0;


      VX /= VW;
      VY /= VW;
      VZ /= VW;

      if (VX < -1.f)
        clip_info |= 1;
      else if (VX > 1.f)
        clip_info |= 2;
      if (VY < -1.f)
        clip_info |= 4;
      else if (VY > 1.f)
        clip_info |= 8;
      if (VZ < -1.f)
        clip_info |= 16;
      else if (VZ > 1.f)
        clip_info |= 32;

      VX = (VX + 1.f)*rd.fb.w * 0.5f;
      VY = (VY + 1.f)*rd.fb.h * 0.5f;
      rd.vertices_x[i] = VX;
      rd.vertices_y[i] = VY;
      rd.vertices_z[i] = VW;
      rd.vertex_clip_info[i] = clip_info;
      }
    }

  template <class TCallBack>
  inline void present(render_data& rd, uint32_t color, TCallBack callback)
    {
    rgb_value rgb;
    rgb.color = color;
    _draw(rd);
    const int border = rd.point_size;

    float inv[16];
    float light[4];
    float tmp[4] = { rd.light_direction[0], rd.light_direction[1], rd.light_direction[2], 0.f };
    invert_orthonormal(inv, rd.camera_position);
    mat_vec_multiply(light, inv, tmp);
#ifdef _AVX2
    __m256 light_avx[3] = { _mm256_set1_ps(light[0]), _mm256_set1_ps(light[1]), _mm256_set1_ps(light[2]) };
#else
    __m128 light_sse[3] = { _mm_set1_ps(light[0]), _mm_set1_ps(light[1]), _mm_set1_ps(light[2]) };
#endif

    float light_color[3];
    light_color[0] = (rd.light_color & 0xff) / 255.f;
    light_color[1] = ((rd.light_color >> 8) & 0xff) / 255.f;
    light_color[2] = ((rd.light_color >> 16) & 0xff) / 255.f;

#ifdef _AVX2
    __m256 light_color_avx[3];
    light_color_avx[0] = _mm256_set1_ps(light_color[0]);
    light_color_avx[1] = _mm256_set1_ps(light_color[1]);
    light_color_avx[2] = _mm256_set1_ps(light_color[2]);
    uint32_t sz = rd.obj.number_of_vertices - (rd.obj.number_of_vertices & 7);
#else
    __m128 light_color_sse[3];
    light_color_sse[0] = _mm_set1_ps(light_color[0]);
    light_color_sse[1] = _mm_set1_ps(light_color[1]);
    light_color_sse[2] = _mm_set1_ps(light_color[2]);
    uint32_t sz = rd.obj.number_of_vertices - (rd.obj.number_of_vertices & 3);
#endif


    const float* p_vert_x = (const float*)rd.vertices_x.data();
    const float* p_vert_y = (const float*)rd.vertices_y.data();
    const float* p_vert_z = (const float*)rd.vertices_z.data();

#ifdef _AVX2
    __m256i colors = _mm256_set1_epi32(color);

    const __m256i zero = _mm256_set1_epi32(0);
    const __m256i wh = _mm256_set1_epi32(rd.fb.w*rd.fb.h);
    const __m256i w_epi32 = _mm256_set1_epi32(rd.fb.w);
    const __m256 halff = _mm256_set1_ps(0.5f);
    const __m256 zerof = _mm256_set1_ps(0.f);
    const __m256 onef = _mm256_set1_ps(1.f);
    const __m256i borderi = _mm256_set1_epi32(border);
    const __m256i w_minus_borderi = _mm256_set1_epi32(rd.fb.w - 1 - border);
    const __m256i h_minus_borderi = _mm256_set1_epi32(rd.fb.h - 1 - border);

    const __m256i ff0 = _mm256_set1_epi32(0x000000ff);
    const __m256i ff1 = _mm256_set1_epi32(0x0000ff00);
    const __m256i ff2 = _mm256_set1_epi32(0x00ff0000);
    const __m256i ff3 = _mm256_set1_epi32(0xff000000);
    for (uint32_t i = 0; i < sz; i += 8)
      {
      const float* px = p_vert_x + i;
      const __m256 *m = (const __m256*)px;
      __m256 x = m[0];
      const float* py = p_vert_y + i;
      m = (const __m256*)py;
      __m256 y = m[0];
      const float* pz = p_vert_z + i;
      m = (const __m256*)pz;
      __m256 z = m[0];


      __m256i X = _mm256_cvtps_epi32(x);
      __m256i Y = _mm256_cvtps_epi32(y);

      __m256i mask = _mm256_cmpgt_epi32(borderi, X);
      mask = _mm256_or_si256(mask, _mm256_cmpgt_epi32(borderi, Y));
      mask = _mm256_or_si256(mask, _mm256_cmpgt_epi32(X, w_minus_borderi));
      mask = _mm256_or_si256(mask, _mm256_cmpgt_epi32(Y, h_minus_borderi));

      if (_mm256_test_all_ones(mask))
        continue;

      if (rd.obj.colors)
        {
        const uint32_t* p = rd.obj.colors + i;
        const __m256i* tmp = (const __m256i*) p;
        colors = tmp[0];
        }

      if (rd.obj.normals)
        {
        const float* p = rd.obj.normals + i * 3;
        const __m128* m = (const __m128*) p;
        __m256 m03 = _mm256_castps128_ps256(m[0]); // load lower halves
        __m256 m14 = _mm256_castps128_ps256(m[1]);
        __m256 m25 = _mm256_castps128_ps256(m[2]);
        m03 = _mm256_insertf128_ps(m03, m[3], 1);  // load upper halves
        m14 = _mm256_insertf128_ps(m14, m[4], 1);
        m25 = _mm256_insertf128_ps(m25, m[5], 1);

        const __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE(2, 1, 3, 2)); // upper x's and y's 
        const __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE(1, 0, 2, 1)); // lower y's and z's
        const __m256 nx = _mm256_shuffle_ps(m03, xy, _MM_SHUFFLE(2, 0, 3, 0));
        const __m256 ny = _mm256_shuffle_ps(yz, xy, _MM_SHUFFLE(3, 1, 2, 0));
        const __m256 nz = _mm256_shuffle_ps(yz, m25, _MM_SHUFFLE(3, 0, 3, 1));


        __m256 diffuse = _mm256_add_ps(halff, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(nx, light_avx[0]), _mm256_mul_ps(ny, light_avx[1])), _mm256_mul_ps(nz, light_avx[2])));
        const __m256 mask1 = _mm256_cmp_ps(diffuse, zerof, 1);
        const __m256 mask2 = _mm256_cmp_ps(onef, diffuse, 1);
        diffuse = _mm256_blendv_ps(diffuse, zerof, mask1);
        diffuse = _mm256_blendv_ps(diffuse, onef, mask2);
        const __m256 int_r = _mm256_mul_ps(diffuse, light_color_avx[0]);
        const __m256 int_g = _mm256_mul_ps(diffuse, light_color_avx[1]);
        const __m256 int_b = _mm256_mul_ps(diffuse, light_color_avx[2]);

        const __m256 red = _mm256_cvtepi32_ps(_mm256_and_si256(colors, ff0));
        const __m256 green = _mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_and_si256(colors, ff1), 8));
        const __m256 blue = _mm256_cvtepi32_ps(_mm256_srai_epi32(_mm256_and_si256(colors, ff2), 16));
        const __m256i red2 = _mm256_min_epi32(_mm256_set1_epi32(255), _mm256_cvtps_epi32(_mm256_mul_ps(int_r, red)));
        const __m256i green2 = _mm256_min_epi32(_mm256_set1_epi32(255), _mm256_cvtps_epi32(_mm256_mul_ps(int_g, green)));
        const __m256i blue2 = _mm256_min_epi32(_mm256_set1_epi32(255), _mm256_cvtps_epi32(_mm256_mul_ps(int_b, blue)));
        colors = _mm256_add_epi32(_mm256_add_epi32(_mm256_add_epi32(ff3, _mm256_slli_epi32(blue2, 16)), _mm256_slli_epi32(green2, 8)), red2);
        }

      __m256 depth = _mm256_div_ps(onef, z);
      for (int xx = -border; xx <= border; ++xx)
        {
        for (int yy = -border; yy <= border; ++yy)
          {
          __m256i index = _mm256_add_epi32(_mm256_add_epi32(X, _mm256_set1_epi32(xx)), _mm256_mullo_epi32(w_epi32, _mm256_add_epi32(Y, _mm256_set1_epi32(yy))));
          index = _mm256_max_epi32(zero, index);
          index = _mm256_min_epi32(wh, index);

          __m256 previous_depth = _mm256_i32gather_ps(rd.fb.zbuffer, index, 4);
          __m256i previous_colors = _mm256_i32gather_epi32((const int*)rd.fb.pixels, index, 4);

          __m256 depth_mask = _mm256_cmp_ps(depth, previous_depth, 13);
          __m256i final_mask = _mm256_andnot_si256(mask, _mm256_castps_si256(depth_mask));

          __m256 current_depth = _mm256_blendv_ps(previous_depth, depth, _mm256_castsi256_ps(final_mask));
          __m256i current_colors = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(previous_colors), _mm256_castsi256_ps(colors), _mm256_castsi256_ps(final_mask)));


          const uint32_t i0 = _mm256_extract_epi32(index, 0);
          const uint32_t i1 = _mm256_extract_epi32(index, 1);
          const uint32_t i2 = _mm256_extract_epi32(index, 2);
          const uint32_t i3 = _mm256_extract_epi32(index, 3);
          const uint32_t i4 = _mm256_extract_epi32(index, 4);
          const uint32_t i5 = _mm256_extract_epi32(index, 5);
          const uint32_t i6 = _mm256_extract_epi32(index, 6);
          const uint32_t i7 = _mm256_extract_epi32(index, 7);

          rd.fb.pixels[i0] = _mm256_extract_epi32(current_colors, 0);
          rd.fb.pixels[i1] = _mm256_extract_epi32(current_colors, 1);
          rd.fb.pixels[i2] = _mm256_extract_epi32(current_colors, 2);
          rd.fb.pixels[i3] = _mm256_extract_epi32(current_colors, 3);
          rd.fb.pixels[i4] = _mm256_extract_epi32(current_colors, 4);
          rd.fb.pixels[i5] = _mm256_extract_epi32(current_colors, 5);
          rd.fb.pixels[i6] = _mm256_extract_epi32(current_colors, 6);
          rd.fb.pixels[i7] = _mm256_extract_epi32(current_colors, 7);

          const __m128 lane0 = _mm256_extractf128_ps(current_depth, 0);
          const __m128 lane1 = _mm256_extractf128_ps(current_depth, 1);
          rd.fb.zbuffer[i0] = _mm_cvtss_f32(lane0);
          rd.fb.zbuffer[i1] = _mm_cvtss_f32(_mm_shuffle_ps(lane0, lane0, 1));
          rd.fb.zbuffer[i2] = _mm_cvtss_f32(_mm_movehl_ps(lane0, lane0));
          rd.fb.zbuffer[i3] = _mm_cvtss_f32(_mm_shuffle_ps(lane0, lane0, 3));
          rd.fb.zbuffer[i4] = _mm_cvtss_f32(lane1);
          rd.fb.zbuffer[i5] = _mm_cvtss_f32(_mm_shuffle_ps(lane1, lane1, 1));
          rd.fb.zbuffer[i6] = _mm_cvtss_f32(_mm_movehl_ps(lane1, lane1));
          rd.fb.zbuffer[i7] = _mm_cvtss_f32(_mm_shuffle_ps(lane1, lane1, 3));

          callback(i, index, final_mask);
          }

        }
      }
#else
    __m128i colors = _mm_set1_epi32(color);

    const __m128i zero = _mm_set1_epi32(0);
    const __m128i wh = _mm_set1_epi32(rd.fb.w*rd.fb.h);
    const __m128i w_epi32 = _mm_set1_epi32(rd.fb.w);
    const __m128 halff = _mm_set1_ps(0.5f);
    const __m128 zerof = _mm_set1_ps(0.f);
    const __m128 onef = _mm_set1_ps(1.f);
    const __m128i borderi = _mm_set1_epi32(border);
    const __m128i w_minus_borderi = _mm_set1_epi32(rd.fb.w - 1 - border);
    const __m128i h_minus_borderi = _mm_set1_epi32(rd.fb.h - 1 - border);

    const __m128i ff0 = _mm_set1_epi32(0x000000ff);
    const __m128i ff1 = _mm_set1_epi32(0x0000ff00);
    const __m128i ff2 = _mm_set1_epi32(0x00ff0000);
    const __m128i ff3 = _mm_set1_epi32(0xff000000);
    for (uint32_t i = 0; i < sz; i += 4)
      {
      const float* px = p_vert_x + i;
      const __m128 *m = (const __m128*)px;
      __m128 x = m[0];
      const float* py = p_vert_y + i;
      m = (const __m128*)py;
      __m128 y = m[0];
      const float* pz = p_vert_z + i;
      m = (const __m128*)pz;
      __m128 z = m[0];

      __m128i X = _mm_cvtps_epi32(x);
      __m128i Y = _mm_cvtps_epi32(y);

      __m128i mask = _mm_cmpgt_epi32(borderi, X);
      mask = _mm_or_si128(mask, _mm_cmpgt_epi32(borderi, Y));
      mask = _mm_or_si128(mask, _mm_cmpgt_epi32(X, w_minus_borderi));
      mask = _mm_or_si128(mask, _mm_cmpgt_epi32(Y, h_minus_borderi));

      if (_mm_test_all_ones(mask))
        continue;

      if (rd.obj.colors)
        {
        const uint32_t* p = rd.obj.colors + i;
        const __m128i* tmp = (const __m128i*) p;
        colors = tmp[0];
        }

      if (rd.obj.normals)
        {
        const float* p = rd.obj.normals + i * 3;
        __m128 x0y0z0x1 = _mm_load_ps(p + 0);
        __m128 y1z1x2y2 = _mm_load_ps(p + 4);
        __m128 z2x3y3z3 = _mm_load_ps(p + 8);
        __m128 x2y2x3y3 = _mm_shuffle_ps(y1z1x2y2, z2x3y3z3, _MM_SHUFFLE(2, 1, 3, 2));
        __m128 y0z0y1z1 = _mm_shuffle_ps(x0y0z0x1, y1z1x2y2, _MM_SHUFFLE(1, 0, 2, 1));
        __m128 nx = _mm_shuffle_ps(x0y0z0x1, x2y2x3y3, _MM_SHUFFLE(2, 0, 3, 0)); // x0x1x2x3
        __m128 ny = _mm_shuffle_ps(y0z0y1z1, x2y2x3y3, _MM_SHUFFLE(3, 1, 2, 0)); // y0y1y2y3
        __m128 nz = _mm_shuffle_ps(y0z0y1z1, z2x3y3z3, _MM_SHUFFLE(3, 0, 3, 1)); // z0z1z2z3        

        __m128 diffuse = _mm_add_ps(halff, _mm_add_ps(_mm_add_ps(_mm_mul_ps(nx, light_sse[0]), _mm_mul_ps(ny, light_sse[1])), _mm_mul_ps(nz, light_sse[2])));
        const __m128 mask1 = _mm_cmplt_ps(diffuse, zerof);
        const __m128 mask2 = _mm_cmplt_ps(onef, diffuse);
        diffuse = _mm_blendv_ps(diffuse, zerof, mask1);
        diffuse = _mm_blendv_ps(diffuse, onef, mask2);
        const __m128 int_r = _mm_mul_ps(diffuse, light_color_sse[0]);
        const __m128 int_g = _mm_mul_ps(diffuse, light_color_sse[1]);
        const __m128 int_b = _mm_mul_ps(diffuse, light_color_sse[2]);

        const __m128 red = _mm_cvtepi32_ps(_mm_and_si128(colors, ff0));
        const __m128 green = _mm_cvtepi32_ps(_mm_srai_epi32(_mm_and_si128(colors, ff1), 8));
        const __m128 blue = _mm_cvtepi32_ps(_mm_srai_epi32(_mm_and_si128(colors, ff2), 16));
        const __m128i red2 = _mm_min_epi32(_mm_set1_epi32(255), _mm_cvtps_epi32(_mm_mul_ps(int_r, red)));
        const __m128i green2 = _mm_min_epi32(_mm_set1_epi32(255), _mm_cvtps_epi32(_mm_mul_ps(int_g, green)));
        const __m128i blue2 = _mm_min_epi32(_mm_set1_epi32(255), _mm_cvtps_epi32(_mm_mul_ps(int_b, blue)));
        colors = _mm_add_epi32(_mm_add_epi32(_mm_add_epi32(ff3, _mm_slli_epi32(blue2, 16)), _mm_slli_epi32(green2, 8)), red2);
        }

      __m128 depth = _mm_div_ps(onef, z);
      for (int xx = -border; xx <= border; ++xx)
        {
        for (int yy = -border; yy <= border; ++yy)
          {
          __m128i index = _mm_add_epi32(_mm_add_epi32(X, _mm_set1_epi32(xx)), _mm_mullo_epi32(w_epi32, _mm_add_epi32(Y, _mm_set1_epi32(yy))));
          index = _mm_max_epi32(zero, index);
          index = _mm_min_epi32(wh, index);

          const uint32_t i0 = _mm_extract_epi32(index, 0);
          const uint32_t i1 = _mm_extract_epi32(index, 1);
          const uint32_t i2 = _mm_extract_epi32(index, 2);
          const uint32_t i3 = _mm_extract_epi32(index, 3);

          __m128 previous_depth = _mm_set_ps(rd.fb.zbuffer[i3], rd.fb.zbuffer[i2], rd.fb.zbuffer[i1], rd.fb.zbuffer[i0]);
          __m128i previous_colors = _mm_set_epi32(rd.fb.pixels[i3], rd.fb.pixels[i2], rd.fb.pixels[i1], rd.fb.pixels[i0]);

          __m128 depth_mask = _mm_cmplt_ps(previous_depth, depth); //ge
          __m128i final_mask = _mm_andnot_si128(mask, _mm_castps_si128(depth_mask));

          __m128 current_depth = _mm_blendv_ps(previous_depth, depth, _mm_castsi128_ps(final_mask));
          __m128i current_colors = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(previous_colors), _mm_castsi128_ps(colors), _mm_castsi128_ps(final_mask)));          

          rd.fb.pixels[i0] = _mm_extract_epi32(current_colors, 0);
          rd.fb.pixels[i1] = _mm_extract_epi32(current_colors, 1);
          rd.fb.pixels[i2] = _mm_extract_epi32(current_colors, 2);
          rd.fb.pixels[i3] = _mm_extract_epi32(current_colors, 3);
         
          rd.fb.zbuffer[i0] = _mm_cvtss_f32(current_depth);
          rd.fb.zbuffer[i1] = _mm_cvtss_f32(_mm_shuffle_ps(current_depth, current_depth, 1));
          rd.fb.zbuffer[i2] = _mm_cvtss_f32(_mm_movehl_ps(current_depth, current_depth));
          rd.fb.zbuffer[i3] = _mm_cvtss_f32(_mm_shuffle_ps(current_depth, current_depth, 3));         

          callback(i, index, final_mask);
          }

        }
      }
#endif

    for (uint32_t i = sz; i < rd.obj.number_of_vertices; ++i)
      {
      if (rd.vertex_clip_info[i])
        continue;
      int X = (int)(rd.vertices_x[i]);
      int Y = (int)(rd.vertices_y[i]);
      if (X < border || Y < border || X > rd.fb.w - 1 - border || Y > rd.fb.h - 1 - border)
        continue;

      if (rd.obj.colors)
        {
        rgb.color = rd.obj.colors[i];
        }
      if (rd.obj.normals)
        {
        float nx = rd.obj.normals[i * 3];
        float ny = rd.obj.normals[i * 3 + 1];
        float nz = rd.obj.normals[i * 3 + 2];

        float diffuse = 0.5f + nx * light[0] + ny * light[1] + nz * light[2];
        diffuse = diffuse < 0.f ? 0.f : (1.f < diffuse) ? 1.f : diffuse; // clamp
        float int_r = diffuse * light_color[0];
        float int_g = diffuse * light_color[1];
        float int_b = diffuse * light_color[2];

        int red = (int)min2(255.f, rgb.rgba.red*int_r);
        int green = (int)min2(255.f, rgb.rgba.green*int_g);
        int blue = (int)min2(255.f, rgb.rgba.blue*int_b);

        rgb.color = 0xff000000 | ((int)blue << 16) | ((int)green << 8) | (int)red;
        }
      float z = 1.f / rd.vertices_z[i];
      for (int xx = -border; xx <= border; ++xx)
        {
        for (int yy = -border; yy <= border; ++yy)
          {
          uint32_t idx = X + xx + rd.fb.w*(Y + yy);
          float previous_depth = rd.fb.zbuffer[idx];
          if (z > previous_depth)
            {
            rd.fb.zbuffer[idx] = z;
            rd.fb.pixels[idx] = rgb.color;
            }
          }
        }
      }
    }


  inline void present(render_data& rd, uint32_t color)
    {
#ifdef _AVX2
    present(rd, color, [](uint32_t, const __m256i&, const __m256i&) {});
#else
    present(rd, color, [](uint32_t, const __m128i&, const __m128i&) {});
#endif
    }

  } // namespace jtk