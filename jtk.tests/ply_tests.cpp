#include "ply_tests.h"
#include "test_assert.h"

#define JTK_PLY_IMPLEMENTATION
#include "../jtk/ply.h"
#include "../jtk/vec.h"

#include <string>
#include <vector>

namespace jtk
  {

  namespace
    {
    int read_vec3_coord(p_ply_argument argument)
      {
      float** p_vertices;
      ply_get_argument_user_data(argument, (void**)&p_vertices, nullptr);
      const auto val = ply_get_argument_value(argument);
      *(*p_vertices) = static_cast<float>(val);
      (*p_vertices) += 3;
      return 1;
      }

    int read_face(p_ply_argument argument)
      {
      uint32_t** p_triangle;
      ply_get_argument_user_data(argument, (void**)&p_triangle, nullptr);
      long length, value_index;
      ply_get_argument_property(argument, nullptr, &length, &value_index);
      if (value_index >= 0 && value_index < 3)
        {
        const auto val = ply_get_argument_value(argument);
        *(*p_triangle) = static_cast<uint32_t>(val);
        ++(*p_triangle);
        }
      return 1;
      }


    static int read_texcoord(p_ply_argument argument)
      {
      float** p_uv;
      ply_get_argument_user_data(argument, (void**)&p_uv, NULL);
      long length, value_index;
      ply_get_argument_property(argument, NULL, &length, &value_index);
      if (value_index >= 0 && value_index < 6)
        {
        double val = ply_get_argument_value(argument);
        *(*p_uv) = (float)val;
        ++(*p_uv);
        }
      if (value_index == (length - 1) && (length != 6))
        {
        for (long j = length; j < 6; ++j)
          {
          *(*p_uv) = (float)0.f;
          ++(*p_uv);
          }
        }
      return 1;
      }

    int c_read_ply(uint32_t* nr_of_vertices,
      float** vertices,
      uint32_t* nr_of_triangles,
      uint32_t** triangles,
      float** uv,
      const char* filename)
      {
      *nr_of_vertices = 0;
      *vertices = nullptr;
      *nr_of_triangles = 0;
      *triangles = nullptr;
      *uv = nullptr;

      const auto ply = ply_open(filename, nullptr, 0, nullptr);
      if (!ply)
        return 0;
      if (!ply_read_header(ply))
        return 0;

      float* p_vertex_pointer_x = nullptr;
      float* p_vertex_pointer_y = nullptr;
      float* p_vertex_pointer_z = nullptr;

      const auto nvertices_x = ply_set_read_cb(ply, "vertex", "x", read_vec3_coord, static_cast<void*>(&p_vertex_pointer_x), 0);
      const auto nvertices_y = ply_set_read_cb(ply, "vertex", "y", read_vec3_coord, static_cast<void*>(&p_vertex_pointer_y), 0);
      const auto nvertices_z = ply_set_read_cb(ply, "vertex", "z", read_vec3_coord, static_cast<void*>(&p_vertex_pointer_z), 0);

      if ((nvertices_x != nvertices_y) || (nvertices_x != nvertices_z))
        {
        ply_close(ply);
        return 0;
        }

      *nr_of_vertices = static_cast<uint32_t>(nvertices_x);

      if (nvertices_x > 0)
        *vertices = static_cast<float*>(malloc(*nr_of_vertices * 3 * sizeof(float)));
      p_vertex_pointer_x = static_cast<float*>(*vertices);
      p_vertex_pointer_y = p_vertex_pointer_x + 1;
      p_vertex_pointer_z = p_vertex_pointer_x + 2;

      uint32_t* p_tria_index = nullptr;

      auto ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", read_face, static_cast<void*>(&p_tria_index), 0);
      if (ntriangles == 0)
        ntriangles = ply_set_read_cb(ply, "face", "vertex_index", read_face, static_cast<void*>(&p_tria_index), 0);

      *nr_of_triangles = static_cast<uint32_t>(ntriangles);

      if (ntriangles > 0)
        *triangles = static_cast<uint32_t*>(malloc(*nr_of_triangles * 3 * sizeof(uint32_t)));

      p_tria_index = static_cast<uint32_t*>(*triangles);

      float* p_uv = nullptr;

      auto ntexcoords = ply_set_read_cb(ply, "face", "texcoord", read_texcoord, (void*)(&p_uv), 0);

      if (ntexcoords > 0)
        *uv = static_cast<float*>(malloc(*nr_of_triangles * 6 * sizeof(float)));

      p_uv = static_cast<float*>(*uv);

      if (!ply_read(ply))
        return 0;

      ply_close(ply);

      return 1;
      }

    int c_read_ply_from_memory(uint32_t* nr_of_vertices,
      float** vertices,
      uint32_t* nr_of_triangles,
      uint32_t** triangles,
      float** uv,
      const char* buffer,
      uint64_t buffer_size)
      {
      *nr_of_vertices = 0;
      *vertices = nullptr;
      *nr_of_triangles = 0;
      *triangles = nullptr;
      *uv = nullptr;

      const auto ply = ply_open_from_memory(buffer, buffer_size, nullptr, 0, nullptr);
      if (!ply)
        return 0;
      if (!ply_read_header(ply))
        return 0;

      float* p_vertex_pointer_x = nullptr;
      float* p_vertex_pointer_y = nullptr;
      float* p_vertex_pointer_z = nullptr;

      const auto nvertices_x = ply_set_read_cb(ply, "vertex", "x", read_vec3_coord, static_cast<void*>(&p_vertex_pointer_x), 0);
      const auto nvertices_y = ply_set_read_cb(ply, "vertex", "y", read_vec3_coord, static_cast<void*>(&p_vertex_pointer_y), 0);
      const auto nvertices_z = ply_set_read_cb(ply, "vertex", "z", read_vec3_coord, static_cast<void*>(&p_vertex_pointer_z), 0);

      if ((nvertices_x != nvertices_y) || (nvertices_x != nvertices_z))
        {
        ply_close(ply);
        return 0;
        }

      *nr_of_vertices = static_cast<uint32_t>(nvertices_x);

      if (nvertices_x > 0)
        *vertices = static_cast<float*>(malloc(*nr_of_vertices * 3 * sizeof(float)));
      p_vertex_pointer_x = static_cast<float*>(*vertices);
      p_vertex_pointer_y = p_vertex_pointer_x + 1;
      p_vertex_pointer_z = p_vertex_pointer_x + 2;

      uint32_t* p_tria_index = nullptr;

      auto ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", read_face, static_cast<void*>(&p_tria_index), 0);
      if (ntriangles == 0)
        ntriangles = ply_set_read_cb(ply, "face", "vertex_index", read_face, static_cast<void*>(&p_tria_index), 0);

      *nr_of_triangles = static_cast<uint32_t>(ntriangles);

      if (ntriangles > 0)
        *triangles = static_cast<uint32_t*>(malloc(*nr_of_triangles * 3 * sizeof(uint32_t)));

      p_tria_index = static_cast<uint32_t*>(*triangles);

      float* p_uv = nullptr;

      auto ntexcoords = ply_set_read_cb(ply, "face", "texcoord", read_texcoord, (void*)(&p_uv), 0);

      if (ntexcoords > 0)
        *uv = static_cast<float*>(malloc(*nr_of_triangles * 6 * sizeof(float)));

      p_uv = static_cast<float*>(*uv);

      if (!ply_read(ply))
        return 0;

      ply_close(ply);

      return 1;
      }

    bool read_ply(std::vector<jtk::vec3<float>>& vertices_vec, std::vector<jtk::vec3<uint32_t>>& triangles_vec, std::vector<jtk::vec3<jtk::vec2<float>>>& uv, const char* filename)
      {
      uint32_t nr_of_vertices;
      float* vertices;
      uint32_t nr_of_triangles;
      uint32_t* triangles;
      float* texcoords;

      if (c_read_ply(&nr_of_vertices, &vertices, &nr_of_triangles, &triangles, &texcoords, filename))
        {
        vertices_vec.resize(nr_of_vertices);
        triangles_vec.resize(nr_of_triangles);
        uv.clear();
        if (texcoords)
          uv.resize(nr_of_triangles);

        for (uint32_t v = 0; v < nr_of_vertices; ++v)
          {
          for (auto j = 0; j < 3; ++j)
            vertices_vec[v][j] = vertices[v * 3 + j];
          }
        for (uint32_t t = 0; t < nr_of_triangles; ++t)
          {
          for (auto j = 0; j < 3; ++j)
            triangles_vec[t][j] = triangles[t * 3 + j];
          if (texcoords)
            {
            for (auto j = 0; j < 3; ++j)
              {
              uv[t][j][0] = texcoords[t * 6 + j * 2];
              uv[t][j][1] = texcoords[t * 6 + j * 2 + 1];
              }
            }
          }

        free(vertices);
        free(triangles);
        free(texcoords);
        return true;
        }
      return false;
      }

    bool read_ply_from_memory(std::vector<jtk::vec3<float>>& vertices_vec, std::vector<jtk::vec3<uint32_t>>& triangles_vec, std::vector<jtk::vec3<jtk::vec2<float>>>& uv, const char* buffer, uint64_t buffer_size)
      {
      uint32_t nr_of_vertices;
      float* vertices;
      uint32_t nr_of_triangles;
      uint32_t* triangles;
      float* texcoords;

      if (c_read_ply_from_memory(&nr_of_vertices, &vertices, &nr_of_triangles, &triangles, &texcoords, buffer, buffer_size))
        {
        vertices_vec.resize(nr_of_vertices);
        triangles_vec.resize(nr_of_triangles);
        uv.clear();
        if (texcoords)
          uv.resize(nr_of_triangles);

        for (uint32_t v = 0; v < nr_of_vertices; ++v)
          {
          for (auto j = 0; j < 3; ++j)
            vertices_vec[v][j] = vertices[v * 3 + j];
          }
        for (uint32_t t = 0; t < nr_of_triangles; ++t)
          {
          for (auto j = 0; j < 3; ++j)
            triangles_vec[t][j] = triangles[t * 3 + j];
          if (texcoords)
            {
            for (auto j = 0; j < 3; ++j)
              {
              uv[t][j][0] = texcoords[t * 6 + j * 2];
              uv[t][j][1] = texcoords[t * 6 + j * 2 + 1];
              }
            }
          }

        free(vertices);
        free(triangles);
        free(texcoords);
        return true;
        }
      return false;
      }


    } // anonymous namespace


  void test_read_ply_file()
    {
    std::string filename = "data/face.ply";
    std::vector<jtk::vec3<float>> vertices;
    std::vector<jtk::vec3<uint32_t>> triangles;
    std::vector<jtk::vec3<jtk::vec2<float>>> uv;
    TEST_ASSERT(read_ply(vertices, triangles, uv, filename.c_str()));
    TEST_EQ((int)1220, (int)vertices.size());
    TEST_EQ((int)2304, (int)triangles.size());
    TEST_EQ((int)0, (int)uv.size());
    }

  void test_read_ply_file_binary()
    {
    std::string filename = "data/face_binary.ply";
    std::vector<jtk::vec3<float>> vertices;
    std::vector<jtk::vec3<uint32_t>> triangles;
    std::vector<jtk::vec3<jtk::vec2<float>>> uv;
    TEST_ASSERT(read_ply(vertices, triangles, uv, filename.c_str()));
    TEST_EQ((int)1220, (int)vertices.size());
    TEST_EQ((int)2304, (int)triangles.size());
    TEST_EQ((int)0, (int)uv.size());
    }

  }

void run_all_ply_tests()
  {
  using namespace jtk;
  test_read_ply_file();
  test_read_ply_file_binary();
  }