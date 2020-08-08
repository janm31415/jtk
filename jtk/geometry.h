#pragma once

#include "concurrency.h"
#include "vec.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace jtk
  {

  /////////////////////////////////////////////////////////////////////////
  // interfaces
  /////////////////////////////////////////////////////////////////////////

  bool read_stl(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  bool write_stl(const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const vec3<float>* triangle_normals, const unsigned short* attributes, const char* filename);
  bool read_obj(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  bool read_obj(std::string& mtl_filename, std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv, const char* filename);
  bool read_texture_filename_from_mtl(std::string& texture_file, const char* filename);
  bool read_off(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  bool write_off(uint32_t nr_of_vertices, const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const char* filename);

  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts);
  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts, const std::vector<vec3<T>>& normals);
  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts, const std::vector<uint32_t>& clrs);
  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts, const std::vector<vec3<T>>& normals, const std::vector<uint32_t>& clrs);
  void write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<uint32_t>>& triangles);
  void write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<uint32_t>& pts_colors, const std::vector<vec3<uint32_t>>& triangles);
    
  class adjacency_list
    {
    public:
      typedef std::vector<uint32_t>::const_iterator const_iterator;

      adjacency_list();
      adjacency_list(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles);
      ~adjacency_list();
      void build(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles);

      const_iterator begin(uint32_t vertex_index) const;
      const_iterator end(uint32_t vertex_index) const;

      uint32_t size() const;
      bool empty() const;

    private:
      std::vector<uint32_t> offset_per_vertex;
      std::vector<uint32_t> triangles_per_vertex_list;
      uint32_t _nr_of_vertices;
    };

  class mutable_adjacency_list
    {
    public:
      typedef std::vector<uint32_t>::const_iterator const_iterator;

      mutable_adjacency_list();
      mutable_adjacency_list(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles);
      ~mutable_adjacency_list();
      void build(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles);
      void add_triangle_to_vertex(uint32_t vertex_id, uint32_t triangle_id);
      void remove_triangle_from_vertex(uint32_t vertex_id, uint32_t triangle_id);

      // triangle_indices should be sorted
      void add_triangles_to_vertex(uint32_t vertex_index, const std::vector<uint32_t>& triangle_indices);
      void remove_triangles_from_vertex(uint32_t vertex_index, const std::vector<uint32_t>& triangle_indices);

      const_iterator begin(uint32_t vertex_index) const;
      const_iterator end(uint32_t vertex_index) const;

      uint32_t size() const;
      uint32_t size(uint32_t vertex_index) const;
      bool empty() const;

    private:
      std::vector<std::vector<uint32_t>> triangles_per_vertex_list;
    };

  template <class TType, class TIndexType>
  void delete_items(std::vector<TType>& vec, const std::vector<TIndexType>& _indices_to_delete);

  void delete_triangles(std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& triangles_to_delete);

  void remove_free_vertices(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles);

  std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  void compute_triangle_normals(std::vector<vec3<float>>& triangle_normals, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles);

  vec3<float> normal(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles);

  vec3<float> normal(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  float area(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles);

  std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list);

  std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list);

  float signed_volume(const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles);

  bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* triangle_normals, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  bool vertex_on_sharp_edge(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  std::vector<std::vector<uint32_t>> triangle_neighbour_indices(const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, const adjacency_list& adj_list);

  // returns false if tria1 and tria2 share three or no edges (in the latter case e contains std::numeric_limits<size_t>::max() twice)
  // otherwise returns the shared edge oriented such as it appears in tria1
  bool edge_between_triangles(uint32_t& v0, uint32_t& v1, const vec3<uint32_t>& tria1, const vec3<uint32_t>& tria2);

  bool is_manifold_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list);

  bool is_boundary_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list);

  bool is_boundary_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list);

  bool is_boundary_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  bool is_boundary_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  enum class shell_connectivity { vertex, edge, manifold };
  std::vector<uint32_t> shells(const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, shell_connectivity conn = shell_connectivity::edge);


  struct trivial_ear
    {
    trivial_ear();
    trivial_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list);
    void compute_quality();
    void compute_angle(const vec3<uint32_t>* triangles);
    bool is_concave() const;
    bool operator < (const trivial_ear& other) const;

    uint32_t v0, v1, v2;
    const vec3<float>* vertices;
    const mutable_adjacency_list* adj_list;
    vec3<float> n;
    float quality;
    float angle_rad;
    };

  inline float dihedral_angle(const vec3<float>& u, const vec3<float>& v, const vec3<float>& a, const vec3<float>& b);

  struct minimum_weight_ear : public trivial_ear
    {
    minimum_weight_ear();
    minimum_weight_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list);
    virtual void compute_quality(const vec3<uint32_t>* triangles);
    virtual inline bool operator < (const minimum_weight_ear& other) const;

    float aspect_ratio;
    float dihedral_rad;
    };

  template <class Ear>
  void fill_hole_ear(std::vector<vec3<uint32_t>>& triangles, const vec3<float>* vertices, mutable_adjacency_list& mal, const std::vector<uint32_t>& hole);
  template <class Ear>
  void fill_holes(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, uint32_t max_holes_size);

  void smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, uint32_t iterations = 100, float lambda = 0.33f, float mu = -0.34f);
  void local_smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& vertex_indices, uint32_t iterations = 100, float lambda = 0.33f, float mu = -0.34f);
  void dyadic_subdivide(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles);

  /////////////////////////////////////////////////////////////////////////
  // implementation
  /////////////////////////////////////////////////////////////////////////



  inline bool read_stl(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
    {
    FILE* inputfile;

    inputfile = fopen(filename, "rb");

    if (!inputfile)
      return false;

    char buffer[80];
    fread(buffer, 1, 80, inputfile);

    if (buffer[0] == 's' && buffer[1] == 'o' && buffer[2] == 'l' && buffer[3] == 'i' && buffer[4] == 'd')
      {
      fclose(inputfile);
      return false;
      }

    uint32_t nr_of_triangles;
    fread((void*)(&nr_of_triangles), sizeof(uint32_t), 1, inputfile);

    {
    std::vector<vec3<float>>().swap(vertices);
    }
    {
    std::vector<vec3<uint32_t>>().swap(triangles);
    }

    std::vector<vec4<float>> vert;
    vert.resize(nr_of_triangles * 3);
    triangles.resize(nr_of_triangles);
    uint32_t count_triangles = 0;
    uint32_t nr_of_vertices = 0;
    auto vert_it = vert.begin();
    auto tria_it = triangles.begin();
    while (!feof(inputfile) && count_triangles < nr_of_triangles)
      {
      ++count_triangles;
      fread(buffer, 1, 50, inputfile);
      (*vert_it)[0] = *reinterpret_cast<float*>(&(buffer[12]));
      (*vert_it)[1] = *reinterpret_cast<float*>(&(buffer[16]));
      (*vert_it)[2] = *reinterpret_cast<float*>(&(buffer[20]));
      (*vert_it++)[3] = *reinterpret_cast<const float*>(&nr_of_vertices);
      (*tria_it)[0] = nr_of_vertices++;
      (*vert_it)[0] = *reinterpret_cast<float*>(&(buffer[24]));
      (*vert_it)[1] = *reinterpret_cast<float*>(&(buffer[28]));
      (*vert_it)[2] = *reinterpret_cast<float*>(&(buffer[32]));
      (*vert_it++)[3] = *reinterpret_cast<const float*>(&nr_of_vertices);
      (*tria_it)[1] = nr_of_vertices++;
      (*vert_it)[0] = *reinterpret_cast<float*>(&(buffer[36]));
      (*vert_it)[1] = *reinterpret_cast<float*>(&(buffer[40]));
      (*vert_it)[2] = *reinterpret_cast<float*>(&(buffer[44]));
      (*vert_it++)[3] = *reinterpret_cast<const float*>(&nr_of_vertices);
      (*tria_it)[2] = nr_of_vertices++;
      ++tria_it;
      }
    fclose(inputfile);
    if (count_triangles != nr_of_triangles)
      return false;

    auto less_fie = [](const vec4<float>& left, const vec4<float>& right)
      {
      return left[0] == right[0] ? (left[1] == right[1] ? left[2] < right[2] : left[1] < right[1]) : left[0] < right[0];
      };

    auto equal_fie = [](const vec4<float>& left, const vec4<float>& right)
      {
      return (left[0] == right[0] && left[1] == right[1] && left[2] == right[2]);
      };

    std::sort(vert.begin(), vert.end(), less_fie);

    nr_of_vertices = 1;
    auto first = vert.begin();
    auto last = vert.end();

    auto result = first;
    while (++first != last)
      {
      if (!equal_fie(*result, *first))
        {
        ++nr_of_vertices;
        result = first;
        }
      }

    vertices.resize(nr_of_vertices);

    nr_of_vertices = 0;
    first = vert.begin();
    result = first;
    vertices[nr_of_vertices][0] = (*result)[0];
    vertices[nr_of_vertices][1] = (*result)[1];
    vertices[nr_of_vertices][2] = (*result)[2];
    triangles[*reinterpret_cast<uint32_t*>(&((*result)[3])) / 3][*reinterpret_cast<uint32_t*>(&((*result)[3])) % 3] = nr_of_vertices;
    while (++first != last)
      {
      if (!equal_fie(*result, *first))
        {
        while (++result != first)
          {
          triangles[*reinterpret_cast<uint32_t*>(&((*result)[3])) / 3][*reinterpret_cast<uint32_t*>(&((*result)[3])) % 3] = nr_of_vertices;
          }
        ++nr_of_vertices;
        vertices[nr_of_vertices][0] = (*result)[0];
        vertices[nr_of_vertices][1] = (*result)[1];
        vertices[nr_of_vertices][2] = (*result)[2];
        triangles[*reinterpret_cast<uint32_t*>(&((*result)[3])) / 3][*reinterpret_cast<uint32_t*>(&((*result)[3])) % 3] = nr_of_vertices;
        }
      }
    while (++result != first)
      {
      triangles[*reinterpret_cast<uint32_t*>(&((*result)[3])) / 3][*reinterpret_cast<uint32_t*>(&((*result)[3])) % 3] = nr_of_vertices;
      }
    assert((++nr_of_vertices) == vertices.size());
    return true;
    }

  inline bool write_stl(const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const vec3<float>* triangle_normals, const unsigned short* attributes, const char* filename)
    {
    FILE* outputfile;
    outputfile = fopen(filename, "wb");


    if (!outputfile)
      return false;

    char buffer[80] = "STL Binary File Format                                                         ";
    fwrite(buffer, 1, 80, outputfile);
    fwrite((void*)(&nr_of_triangles), sizeof(uint32_t), 1, outputfile);
    uint32_t count_triangles = 0;
    const vec3<uint32_t>* tria_it = triangles;
    const vec3<float>* tria_norm_it = triangle_normals;
    const unsigned short* attr_it = attributes;
    while (count_triangles < nr_of_triangles)
      {
      ++count_triangles;
      vec3<uint32_t> tria = *tria_it++;
      if (tria_norm_it)
        {
        vec3<float> tria_norm = *tria_norm_it++;
        memcpy((void*)&buffer[0], &tria_norm[0], sizeof(float));
        memcpy((void*)&buffer[4], &tria_norm[1], sizeof(float));
        memcpy((void*)&buffer[8], &tria_norm[2], sizeof(float));
        }
      else
        {
        float z = 0.f;
        memcpy((void*)&buffer[0], &z, sizeof(float));
        memcpy((void*)&buffer[4], &z, sizeof(float));
        memcpy((void*)&buffer[8], &z, sizeof(float));
        }
      memcpy((void*)&buffer[12], &vertices[tria[0]][0], sizeof(float));
      memcpy((void*)&buffer[16], &vertices[tria[0]][1], sizeof(float));
      memcpy((void*)&buffer[20], &vertices[tria[0]][2], sizeof(float));
      memcpy((void*)&buffer[24], &vertices[tria[1]][0], sizeof(float));
      memcpy((void*)&buffer[28], &vertices[tria[1]][1], sizeof(float));
      memcpy((void*)&buffer[32], &vertices[tria[1]][2], sizeof(float));
      memcpy((void*)&buffer[36], &vertices[tria[2]][0], sizeof(float));
      memcpy((void*)&buffer[40], &vertices[tria[2]][1], sizeof(float));
      memcpy((void*)&buffer[44], &vertices[tria[2]][2], sizeof(float));
      if (!attr_it)
        {
        unsigned short z = 0;
        memcpy((void*)&buffer[48], &z, sizeof(unsigned short));
        }
      else
        {
        ++attr_it;
        memcpy((void*)&buffer[48], attr_it, sizeof(unsigned short));
        }
      fwrite(buffer, 1, 50, outputfile);
      }
    fclose(outputfile);
    return true;
    }

  inline bool read_off(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
    {
    uint32_t nr_of_vertices = 0;
    uint32_t nr_of_triangles = 0;
    vertices.clear();
    triangles.clear();
    FILE* inputfile = fopen(filename, "r");
    if (!inputfile)
      return false;
    char str[80];
    fscanf(inputfile, "%s\n", str);
    if (std::string(str) != "OFF" && std::string(str) != "COFF")
      {
      fclose(inputfile);
      return false;
      }
    bool color = std::string(str) == "COFF";
    uint32_t numedges = 0;
    fscanf(inputfile, "%d %d %d\n", &nr_of_vertices, &nr_of_triangles, &numedges);
    vertices.resize(nr_of_vertices);
    triangles.resize(nr_of_triangles);
    for (auto& v : vertices)
      {
      float x = 0.0;
      float y = 0.0;
      float z = 0.0;
      int r = 0;
      int g = 0;
      int b = 0;
      int a = 0;
      if (color)
        fscanf(inputfile, "%f %f %f %d %d %d %d\n", &x, &y, &z, &r, &g, &b, &a);
      else
        fscanf(inputfile, "%f %f %f\n", &x, &y, &z);
      v[0] = x;
      v[1] = y;
      v[2] = z;
      }
    for (auto& tria : triangles)
      {
      uint32_t polysize = 0;
      uint32_t t0, t1, t2;
      fscanf(inputfile, "%d %d %d %d\n", &polysize, &t0, &t1, &t2);
      if (polysize != 3)
        {
        vertices.clear();
        triangles.clear();
        nr_of_vertices = 0;
        nr_of_triangles = 0;
        fclose(inputfile);
        return false;
        }
      tria[0] = t0;
      tria[1] = t1;
      tria[2] = t2;
      }
    fclose(inputfile);
    return true;
    }

  inline bool write_off(uint32_t nr_of_vertices, const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const char* filename)
    {
    FILE* outputfile = fopen(filename, "w");
    if (!outputfile)
      return false;
    fprintf(outputfile, "OFF\n");
    fprintf(outputfile, "%d %d 0\n", nr_of_vertices, nr_of_triangles);
    // vertex data
    const vec3<float>* p_vert = vertices;
    for (uint32_t i = 0; i < nr_of_vertices; ++i)
      {
      float x = (*p_vert)[0];
      float y = (*p_vert)[1];
      float z = (*p_vert)[2];
      ++p_vert;
      fprintf(outputfile, "%f %f %f\n", x, y, z);
      }
    // face data
    const vec3<uint32_t>* p_tria = triangles;
    for (uint32_t i = 0; i < nr_of_triangles; ++i)
      {
      uint32_t t0 = (*p_tria)[0];
      uint32_t t1 = (*p_tria)[1];
      uint32_t t2 = (*p_tria)[2];
      ++p_tria;
      fprintf(outputfile, "3 %d %d %d\n", t0, t1, t2);
      }
    fclose(outputfile);
    return true;
    }

  inline bool read_obj(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
    {
    FILE* f = nullptr;
    auto ferr = fopen_s(&f, filename, "r");
    if (ferr != 0)
      return false;
    if (!f)
      return false;


    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      if (buffer[0] == 'v' && buffer[1] == ' ')
        {
        float x, y, z;
        auto err = sscanf(buffer, "v %f %f %f\n", &x, &y, &z);
        if (err != 3)
          {
          fclose(f);
          return false;
          }
        vertices.push_back(vec3<float>(x, y, z));
        }
      else if (buffer[0] == 'f' && buffer[1] == ' ')
        {
        uint32_t t0, t1, t2, v0, v1, v2;
        auto err = sscanf(buffer, "f %d/%d %d/%d %d/%d\n", &t0, &v0, &t1, &v1, &t2, &v2);
        if (err == 6)
          {
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          }
        else
          {
          err = sscanf(buffer, "f %d//%d %d//%d %d//%d\n", &t0, &v0, &t1, &v1, &t2, &v2);
          if (err != 6)
            {
            err = sscanf(buffer, "f %d %d %d\n", &t0, &t1, &t2);
            if (err != 3)
              {
              uint32_t tx0, tx1, tx2;
              err = sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", &t0, &v0, &tx0, &t1, &v1, &tx1, &t2, &v2, &tx2);
              if (err != 9)
                {
                fclose(f);
                return false;
                }
              }
            }
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          }
        }
      }
    fclose(f);
    return true;
    }

  inline bool read_texture_filename_from_mtl(std::string& texture_file, const char* filename)
    {
    FILE* f = nullptr;
    auto err = fopen_s(&f, filename, "r");
    if (err != 0)
      return false;
    if (!f)
      return false;
    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      if (buffer[0] == 'm' && buffer[1] == 'a' && buffer[2] == 'p' && buffer[3] == '_' && buffer[4] == 'K' && buffer[5] == 'd')
        {
        char texture[256];
        auto scan_err = sscanf(buffer, "map_Kd %s\n", texture);
        scan_err;
        texture_file = std::string(texture);
        fclose(f);
        return true;
        }
      }
    fclose(f);
    return false;
    }

  inline bool read_obj(std::string& mtl_filename, std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv, const char* filename)
    {
    mtl_filename = "";
    std::vector<vec3<float>>().swap(vertices);
    std::vector<vec3<uint32_t>>().swap(triangles);
    std::vector<vec3<vec2<float>>>().swap(uv);
    FILE* f = nullptr;
    auto ferr = fopen_s(&f, filename, "r");
    if (ferr != 0)
      return false;
    if (!f)
      return false;

    std::vector<vec2<float>> tex;
    std::vector<vec3<uint32_t>> tria_uv;

    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      if (buffer[0] == 'm' && buffer[1] == 't' && buffer[2] == 'l' && buffer[3] == 'l' && buffer[4] == 'i' && buffer[5] == 'b')
        {
        char mtlfile[256];
        auto scan_err = sscanf(buffer, "mtllib %s\n", mtlfile);
        if (scan_err == 1)
          {
          mtl_filename = std::string(mtlfile);
          }
        }
      if (buffer[0] == 'v' && buffer[1] == ' ')
        {
        float x, y, z;
        auto err = sscanf(buffer, "v %f %f %f\n", &x, &y, &z);
        if (err != 3)
          {
          fclose(f);
          return false;
          }
        vertices.push_back(vec3<float>(x, y, z));
        }
      else if (buffer[0] == 'v' && buffer[1] == 't' && buffer[2] == ' ')
        {
        float x, y;
        auto err = sscanf(buffer, "vt %f %f\n", &x, &y);
        if (err != 2)
          {
          fclose(f);
          return false;
          }
        tex.push_back(vec2<float>(x, y));
        }
      else if (buffer[0] == 'f' && buffer[1] == ' ')
        {
        uint32_t t0, t1, t2, v0, v1, v2;
        auto err = sscanf(buffer, "f %d/%d %d/%d %d/%d\n", &t0, &v0, &t1, &v1, &t2, &v2);
        if (err == 6)
          {
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          }
        else
          {
          err = sscanf(buffer, "f %d//%d %d//%d %d//%d\n", &t0, &v0, &t1, &v1, &t2, &v2);
          if (err == 6)
            {
            tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            }
          else
            {
            err = sscanf(buffer, "f %d %d %d\n", &t0, &t1, &t2);
            if (err == 3)
              {
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              }
            else
              {
              uint32_t tx0, tx1, tx2;
              err = sscanf(buffer, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", &t0, &v0, &tx0, &t1, &v1, &tx1, &t2, &v2, &tx2);
              if (err != 9)
                {
                fclose(f);
                return false;
                }
              tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              }
            }
          }
        }
      }
    fclose(f);
    if (!tria_uv.empty() && (triangles.size() != tria_uv.size()))
      return false;
    if (!tria_uv.empty())
      {
      for (auto& t : tex)
        {
        t[1] = 1.f - t[1];
        }
      uv.reserve(triangles.size());
      for (auto t : tria_uv)
        {
        uv.emplace_back(tex[t[0]], tex[t[1]], tex[t[2]]);
        }
      }
    return true;
    }

  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts)
    {
    FILE* fp;
    fopen_s(&fp, filename, "wt");
    fprintf(fp, "ply\n");
    fprintf(fp, "format ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "end_header\n");
    for (int i = 0; i < pts.size(); ++i)
      {
      auto point = pts[i];
      fprintf(fp, "%f %f %f\n", (float)point[0], (float)point[1], (float)point[2]);
      }
    fclose(fp);
    }

  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts, const std::vector<vec3<T>>& normals)
    {
    if (normals.empty())
      {
      write_ply(filename, pts);
      return;
      }
    FILE* fp;
    fopen_s(&fp, filename, "wt");
    fprintf(fp, "ply\n");
    fprintf(fp, "format ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "property float nx\n");
    fprintf(fp, "property float ny\n");
    fprintf(fp, "property float nz\n");
    fprintf(fp, "end_header\n");
    for (int i = 0; i < pts.size(); ++i)
      {
      auto point = pts[i];
      auto norm = normals[i];
      fprintf(fp, "%f %f %f %f %f %f\n", (float)point[0], (float)point[1], (float)point[2], (float)norm[0], (float)norm[1], (float)norm[2]);
      }
    fclose(fp);
    }  

  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts, const std::vector<uint32_t>& clrs)
    {
    if (clrs.empty())
      {
      write_ply(filename, pts);
      return;
      }
    FILE* fp;
    fopen_s(&fp, filename, "wt");
    fprintf(fp, "ply\n");
    fprintf(fp, "format ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "property uchar red\n");
    fprintf(fp, "property uchar green\n");
    fprintf(fp, "property uchar blue\n");
    fprintf(fp, "end_header\n");
    for (int i = 0; i < pts.size(); ++i)
      {
      auto point = pts[i];
      auto clr = clrs[i];
      fprintf(fp, "%f %f %f %d %d %d\n", (float)point[0], (float)point[1], (float)point[2], clr & 0xff, (clr & 0xff00) >> 8, (clr & 0xff0000) >> 16);
      }
    fclose(fp);
    }

  template <class T>
  void write_ply(const char* filename, const std::vector<vec3<T>>& pts, const std::vector<vec3<T>>& normals, const std::vector<uint32_t>& clrs)
    {
    if (normals.empty())
      {
      write_ply(filename, pts, clrs);
      return;
      }
    if (clrs.empty())
      {
      write_ply(filename, pts, normals);
      return;
      }
    FILE* fp;
    fopen_s(&fp, filename, "wt");
    fprintf(fp, "ply\n");
    fprintf(fp, "format ascii 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "property float nx\n");
    fprintf(fp, "property float ny\n");
    fprintf(fp, "property float nz\n");
    fprintf(fp, "property uchar red\n");
    fprintf(fp, "property uchar green\n");
    fprintf(fp, "property uchar blue\n");
    fprintf(fp, "end_header\n");
    for (int i = 0; i < pts.size(); ++i)
      {
      auto point = pts[i];
      auto norm = normals[i];
      auto clr = clrs[i];
      fprintf(fp, "%f %f %f %f %f %f %d %d %d\n", (float)point[0], (float)point[1], (float)point[2], (float)norm[0], (float)norm[1], (float)norm[2], clr & 0xff, (clr & 0xff00) >> 8, (clr & 0xff0000) >> 16);
      }
    fclose(fp);
    }

  inline void write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<uint32_t>>& triangles)
    {
    FILE* fp;
    fopen_s(&fp, filename, "wb");
    fprintf(fp, "ply\n");
    fprintf(fp, "format binary_little_endian 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "element face %d\n", (int)triangles.size());
    fprintf(fp, "property list uchar int vertex_indices\n");
    fprintf(fp, "end_header\n");
    fwrite(pts.data(), sizeof(vec3<float>), pts.size(), fp);
    unsigned char tria_size = 3;
    for (int i = 0; i < triangles.size(); ++i)
      {
      fwrite(&tria_size, 1, 1, fp);
      fwrite(triangles.data() + i, sizeof(vec3<uint32_t>), 1, fp);
      }
    fclose(fp);
    }


  inline void write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<uint32_t>& pts_colors, const std::vector<vec3<uint32_t>>& triangles)
    {
    FILE* fp;
    fopen_s(&fp, filename, "wb");
    fprintf(fp, "ply\n");
    fprintf(fp, "format binary_little_endian 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "property uchar red\n");
    fprintf(fp, "property uchar green\n");
    fprintf(fp, "property uchar blue\n");
    fprintf(fp, "element face %d\n", (int)triangles.size());
    fprintf(fp, "property list uchar int vertex_indices\n");
    fprintf(fp, "end_header\n");
    for (int i = 0; i < pts.size(); ++i)
      {
      fwrite(pts.data() + i, sizeof(vec3<float>), 1, fp);
      uint32_t clr = pts_colors[i];
      unsigned char* rgb = reinterpret_cast<unsigned char*>(&clr);
      fwrite(rgb, 1, 3, fp);
      }
    unsigned char tria_size = 3;
    for (int i = 0; i < triangles.size(); ++i)
      {
      fwrite(&tria_size, 1, 1, fp);
      fwrite(triangles.data() + i, sizeof(vec3<uint32_t>), 1, fp);
      }
    fclose(fp);
    }

  inline adjacency_list::adjacency_list() : _nr_of_vertices(0)
    {

    }

  inline adjacency_list::adjacency_list(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles) : _nr_of_vertices(nr_of_vertices)
    {
    build(nr_of_vertices, triangles, nr_of_triangles);
    }

  inline adjacency_list::~adjacency_list()
    {

    }

  inline void adjacency_list::build(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
    {
    _nr_of_vertices = nr_of_vertices;
    std::vector<uint32_t>(nr_of_vertices + 1, 0).swap(offset_per_vertex);
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      ++offset_per_vertex[triangles[t][0]];
      ++offset_per_vertex[triangles[t][1]];
      ++offset_per_vertex[triangles[t][2]];
      }
    uint32_t prev_value = offset_per_vertex[0];
    offset_per_vertex[0] = 0;
    for (uint32_t v = 0; v < nr_of_vertices; ++v)
      {
      uint32_t current = offset_per_vertex[v + 1];
      offset_per_vertex[v + 1] = offset_per_vertex[v] + prev_value;
      prev_value = current;
      }
    std::vector<uint32_t>(offset_per_vertex.back()).swap(triangles_per_vertex_list);
    std::vector<uint32_t> vertex_count(nr_of_vertices, 0);
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      for (int j = 0; j < 3; ++j)
        {
        const uint32_t v = triangles[t][j];
        const uint32_t pos = offset_per_vertex[v];
        const uint32_t off = vertex_count[v]++;
        triangles_per_vertex_list[pos + off] = t;
        }
      }
    }

  inline adjacency_list::const_iterator adjacency_list::begin(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list.begin() + offset_per_vertex[vertex_index];
    }

  inline adjacency_list::const_iterator adjacency_list::end(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list.begin() + offset_per_vertex[vertex_index + 1];
    }

  inline uint32_t adjacency_list::size() const
    {
    return _nr_of_vertices;
    }

  inline bool adjacency_list::empty() const
    {
    return _nr_of_vertices == 0;
    }

  inline mutable_adjacency_list::mutable_adjacency_list()
    {

    }

  inline mutable_adjacency_list::mutable_adjacency_list(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
    {
    build(nr_of_vertices, triangles, nr_of_triangles);
    }

  inline mutable_adjacency_list::~mutable_adjacency_list()
    {

    }

  inline void mutable_adjacency_list::build(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
    {
    triangles_per_vertex_list.resize(nr_of_vertices);
    for (uint32_t v = 0; v < nr_of_vertices; ++v)
      {
      triangles_per_vertex_list[v].reserve(6);
      }
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      triangles_per_vertex_list[triangles[t][0]].push_back(t);
      triangles_per_vertex_list[triangles[t][1]].push_back(t);
      triangles_per_vertex_list[triangles[t][2]].push_back(t);
      }
    }

  inline void mutable_adjacency_list::add_triangle_to_vertex(uint32_t vertex_id, uint32_t triangle_id)
    {
    auto first = triangles_per_vertex_list[vertex_id].begin();
    auto last = triangles_per_vertex_list[vertex_id].end();
    auto it = std::lower_bound(first, last, triangle_id);
    triangles_per_vertex_list[vertex_id].insert(it, triangle_id);
    }

  inline void mutable_adjacency_list::remove_triangle_from_vertex(uint32_t vertex_id, uint32_t triangle_id)
    {
    auto first = triangles_per_vertex_list[vertex_id].begin();
    auto last = triangles_per_vertex_list[vertex_id].end();
    auto it = std::lower_bound(first, last, triangle_id);
    if (it != last && *it == triangle_id)
      {
      triangles_per_vertex_list[vertex_id].erase(it);
      }
    }


  inline void mutable_adjacency_list::add_triangles_to_vertex(uint32_t vertex_index, const std::vector<uint32_t>& triangle_indices)
    {
    auto sz = size(vertex_index);
    triangles_per_vertex_list[vertex_index].resize(sz + triangle_indices.size());
    auto first = triangles_per_vertex_list[vertex_index].begin();
    auto last = first + sz;
    size_t offset = triangle_indices.size();
    auto first_ind = triangle_indices.rbegin();
    while (offset && (first != last))
      {
      if (*(last - 1) < *first_ind)
        {
        *(last + offset - 1) = *first_ind;
        ++first_ind;
        --offset;
        }
      else
        {
        *(last + offset - 1) = *(last - 1);
        --last;
        }
      }
    while (offset)
      {
      *(last + offset - 1) = *first_ind++;
      --offset;
      }
    }

  inline void mutable_adjacency_list::remove_triangles_from_vertex(uint32_t vertex_index, const std::vector<uint32_t>& triangle_indices)
    {
    auto first = triangles_per_vertex_list[vertex_index].begin();
    auto last = triangles_per_vertex_list[vertex_index].end();

    size_t offset = 0;

    auto first_ind = triangle_indices.begin();
    auto last_ind = triangle_indices.end();
    for (; (first + offset) < last; ++first)
      {
      while ((first_ind != last_ind) && (*(first + offset) > *first_ind))
        ++first_ind;
      while ((first_ind != last_ind) && (*(first + offset) == *first_ind))
        {
        ++offset;
        ++first_ind;
        }
      if ((first + offset) < last)
        *first = *(first + offset);
      else
        break;
      }
    triangles_per_vertex_list[vertex_index].resize(triangles_per_vertex_list[vertex_index].size() - offset);
    }

  inline mutable_adjacency_list::const_iterator mutable_adjacency_list::begin(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list[vertex_index].begin();
    }

  inline mutable_adjacency_list::const_iterator mutable_adjacency_list::end(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list[vertex_index].end();
    }

  inline uint32_t mutable_adjacency_list::size(uint32_t vertex_index) const
    {
    return (uint32_t)triangles_per_vertex_list[vertex_index].size();
    }

  inline uint32_t mutable_adjacency_list::size() const
    {
    return (uint32_t)triangles_per_vertex_list.size();
    }

  inline bool mutable_adjacency_list::empty() const
    {
    return triangles_per_vertex_list.empty();
    }


  template <class TType, class TIndexType>
  void delete_items(std::vector<TType>& vec, const std::vector<TIndexType>& _indices_to_delete)
    {
    if (_indices_to_delete.empty() || vec.empty())
      return;
    std::vector<TIndexType> indices_to_delete(_indices_to_delete);
    assert(vec.size() > *std::max_element(indices_to_delete.begin(), indices_to_delete.end()));
    std::sort(indices_to_delete.begin(), indices_to_delete.end());
    indices_to_delete.erase(std::unique(indices_to_delete.begin(), indices_to_delete.end()), indices_to_delete.end());

    if (indices_to_delete.size() == vec.size())
      {
      vec.clear();
      return;
      }

    auto last = --vec.end();

    for (auto rit = indices_to_delete.rbegin(); rit != indices_to_delete.rend(); ++rit)
      {
      TIndexType index = *rit;
      auto it = vec.begin() + index;
      if (it != last)
        {
        std::swap(*it, *last);
        }
      --last;
      }
    vec.erase(++last, vec.end());
    }

  inline void delete_triangles(std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& triangles_to_delete)
    {
    delete_items(triangles, triangles_to_delete);
    }

  inline void remove_free_vertices(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles)
    {
    std::vector<uint32_t> count_occurence(vertices.size(), 0);
    for (const auto& tria : triangles)
      {
      ++count_occurence[tria[0]];
      ++count_occurence[tria[1]];
      ++count_occurence[tria[2]];
      }
    std::vector<uint32_t> diff(vertices.size(), 0);
    uint32_t current_off = 0;
    std::vector<vec3<float>> new_vertices;
    auto cnt = (size_t)std::count(count_occurence.begin(), count_occurence.end(), (uint32_t)0);
    if (cnt <= vertices.size())
      {
      new_vertices.reserve((size_t)(vertices.size() - cnt));
      }
    for (size_t v = 0; v < vertices.size(); ++v)
      {
      if (count_occurence[v] == 0)
        {
        ++current_off;
        }
      else
        new_vertices.push_back(vertices[v]);
      diff[v] = current_off;
      }

    for (auto& tria : triangles)
      {
      tria[0] -= diff[tria[0]];
      tria[1] -= diff[tria[1]];
      tria[2] -= diff[tria[2]];
      }
    vertices.swap(new_vertices);
    }

  inline std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    std::vector<uint32_t> one_ring;
    one_ring.reserve(16);
    auto it = adj_list.begin(vertex_index);
    const auto it_end = adj_list.end(vertex_index);
    for (; it != it_end; ++it)
      {
      const auto tria = triangles[*it];
      uint32_t v = 0;
      if (tria[v] != vertex_index)
        {
        ++v;
        if (tria[v] != vertex_index)
          ++v;
        }
      one_ring.push_back(tria[(v + 1) % 3]);
      one_ring.push_back(tria[(v + 2) % 3]);
      }
    std::sort(one_ring.begin(), one_ring.end());
    one_ring.erase(std::unique(one_ring.begin(), one_ring.end()), one_ring.end());
    return one_ring;
    }

  inline std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    std::vector<uint32_t> one_ring;
    one_ring.reserve(16);
    auto it = adj_list.begin(vertex_index);
    const auto it_end = adj_list.end(vertex_index);
    for (; it != it_end; ++it)
      {
      const auto tria = triangles[*it];
      uint32_t v = 0;
      if (tria[v] != vertex_index)
        {
        ++v;
        if (tria[v] != vertex_index)
          ++v;
        }
      one_ring.push_back(tria[(v + 1) % 3]);
      one_ring.push_back(tria[(v + 2) % 3]);
      }
    std::sort(one_ring.begin(), one_ring.end());
    one_ring.erase(std::unique(one_ring.begin(), one_ring.end()), one_ring.end());
    return one_ring;
    }

  inline vec3<float> normal(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles)
    {
    const auto v0 = triangles[triangle_index][0];
    const auto v1 = triangles[triangle_index][1];
    const auto v2 = triangles[triangle_index][2];
    const vec3<float> V0(vertices[v0]);
    const vec3<float> V1(vertices[v1]);
    const vec3<float> V2(vertices[v2]);
    const auto n = normalize(cross(V1 - V0, V2 - V0));
    return n;
    }

  inline vec3<float> normal(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle)
    {
    vec3<float> n(0, 0, 0);
    auto indices_first = adj_list.begin(vertex_index);
    auto indices_last = adj_list.end(vertex_index);
    if (indices_first == indices_last)
      return n;
    auto ref_n = normal(*indices_first, vertices, triangles);
    uint32_t num_trias = 0;
    for (; indices_first != indices_last; ++indices_first)
      {
      auto tria_n = normal(*indices_first, vertices, triangles);
      if (dot(tria_n, ref_n) >= cos_sharp_angle)
        {
        ++num_trias;
        n[0] += tria_n[0];
        n[1] += tria_n[1];
        n[2] += tria_n[2];
        }
      }
    n[0] /= (float)num_trias;
    n[1] /= (float)num_trias;
    n[2] /= (float)num_trias;
    return n;
    }

  inline void compute_triangle_normals(std::vector<vec3<float>>& triangle_normals, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles)
    {
        {
        std::vector<vec3<float>>().swap(triangle_normals);
        }
        triangle_normals.resize(nr_of_triangles);
        parallel_for(uint32_t(0), nr_of_triangles, [&](uint32_t t)
          {
          const auto v0 = triangles[t][0];
          const auto v1 = triangles[t][1];
          const auto v2 = triangles[t][2];
          const vec3<float> V0(vertices[v0]);
          const vec3<float> V1(vertices[v1]);
          const vec3<float> V2(vertices[v2]);
          const auto n = normalize(cross(V1 - V0, V2 - V0));
          triangle_normals[t] = n;
          });
    }

  inline std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list)
    {
    assert(v0 != v1);
    std::vector<uint32_t> result;
    std::set_intersection(adj_list.begin(v0), adj_list.end(v0), adj_list.begin(v1), adj_list.end(v1), std::back_inserter(result));
    return result;
    }

  inline std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list)
    {
    assert(v0 != v1);
    std::vector<uint32_t> result;
    result.reserve(2);
    std::set_intersection(adj_list.begin(v0), adj_list.end(v0), adj_list.begin(v1), adj_list.end(v1), std::back_inserter(result));
    return result;
    }

  inline float signed_volume(const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles)
    {
    float vol = 0.f;
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      const auto& a = vertices[triangles[t][0]];
      const auto& b = vertices[triangles[t][1]];
      const auto& c = vertices[triangles[t][2]];
      vol += dot(a - vertices[0], cross(b - vertices[0], c - vertices[0])) / 6;
      }
    return vol;
    }

  inline float area(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles)
    {
    const auto v0 = triangles[triangle_index][0];
    const auto v1 = triangles[triangle_index][1];
    const auto v2 = triangles[triangle_index][2];
    const vec3<float> V0(vertices[v0]);
    const vec3<float> V1(vertices[v1]);
    const vec3<float> V2(vertices[v2]);
    return length(cross(V1 - V0, V2 - V0)) / 2.f;
    }

  inline bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle)
    {
    auto trias = triangle_indices_from_edge(v0, v1, adj_list);
    if (trias.size() == 2)
      {
      vec3<float> n[2];
      for (int i = 0; i < 2; ++i)
        {
        const auto vv0 = triangles[trias[i]][0];
        const auto vv1 = triangles[trias[i]][1];
        const auto vv2 = triangles[trias[i]][2];
        const vec3<float> V0(vertices[vv0]);
        const vec3<float> V1(vertices[vv1]);
        const vec3<float> V2(vertices[vv2]);
        n[i] = normalize(cross(V1 - V0, V2 - V0));
        }
      auto ca = dot(n[0], n[1]);
      return ca < cos_sharp_angle;
      }
    return false;
    }

  inline bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* triangle_normals, const adjacency_list& adj_list, float cos_sharp_angle)
    {
    auto trias = triangle_indices_from_edge(v0, v1, adj_list);
    if (trias.size() == 2)
      {
      auto ca = dot(triangle_normals[trias[0]], triangle_normals[trias[1]]);
      return ca < cos_sharp_angle;
      }
    return false;
    }

  inline bool vertex_on_sharp_edge(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle)
    {
    auto one_ring = one_ring_vertices_from_vertex(vertex_index, adj_list, triangles);
    for (auto v : one_ring)
      {
      if (is_sharp_edge(vertex_index, v, vertices, triangles, adj_list, cos_sharp_angle))
        return true;
      }
    return false;
    }

  inline std::vector<std::vector<uint32_t>> triangle_neighbour_indices(const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, const adjacency_list& adj_list)
    {
    std::vector<std::vector<uint32_t>> result(nr_of_triangles);
    parallel_for(uint32_t(0), nr_of_triangles, [&](uint32_t idx)
      {
      auto tria = triangles[idx];
      std::set<uint32_t> triangle_neighbours;
      for (uint32_t i = 0; i < 3; ++i)
        {
        auto edge_tria = triangle_indices_from_edge(tria[(i + 1) % 3], tria[(i + 2) % 3], adj_list);
        for (auto t : edge_tria)
          {
          if (t != idx)
            triangle_neighbours.insert(t);
          }
        }
      result[idx] = std::vector<uint32_t>(triangle_neighbours.begin(), triangle_neighbours.end());
      });
    return result;
    }

  inline bool edge_between_triangles(uint32_t& v0, uint32_t& v1, const vec3<uint32_t>& tria1, const vec3<uint32_t>& tria2)
    {
    uint32_t idx = 0;
    std::array<uint32_t, 3> common, idces;
    std::array<bool, 3> found = { false, false, false };
    for (uint32_t i = 0; i < 3; ++i)
      for (uint32_t j = 0; j < 3; ++j)
        if (!found[j] && tria1[i] == tria2[j])
          {
          common[idx] = tria1[i];
          idces[idx] = i;
          found[j] = true;
          ++idx;
          break;
          }
    if (idx < 2)
      {
      v0 = v1 = std::numeric_limits<uint32_t>::max();
      return false;
      }
    v0 = common[0];
    v1 = common[1];
    if ((idces[0] + 1) % 3 != idces[1])
      std::swap(v0, v1);
    return (idx == 2);
    }

  inline bool is_manifold_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list)
    {
    auto tri_idces = triangle_indices_from_edge(v0, v1, adj_list);
    return tri_idces.size() == 2;
    }

  inline bool is_boundary_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list)
    {
    return triangle_indices_from_edge(v0, v1, adj_list).size() == 1;
    }

  inline bool is_boundary_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list)
    {
    return triangle_indices_from_edge(v0, v1, adj_list).size() == 1;
    }

  inline bool is_boundary_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    auto vert = one_ring_vertices_from_vertex(vertex_index, adj_list, triangles);
    for (auto v : vert)
      {
      if (is_boundary_edge(v, vertex_index, adj_list))
        return true;
      }
    return false;
    }

  inline bool is_boundary_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    auto vert = one_ring_vertices_from_vertex(vertex_index, adj_list, triangles);
    for (auto v : vert)
      {
      if (is_boundary_edge(v, vertex_index, adj_list))
        return true;
      }
    return false;
    }


  namespace details
    {
    inline std::vector<uint32_t> _shells(const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, bool manifold)
      {
      std::vector<uint32_t> sh(nr_of_triangles, 0);
      auto adj_list = adjacency_list(nr_of_vertices, triangles, nr_of_triangles);
      auto tr_neigh = triangle_neighbour_indices(triangles, nr_of_triangles, adj_list);
      std::vector<bool> marked(nr_of_triangles, false);
      uint32_t current_shell = 0;
      for (uint32_t tria_index = 0; tria_index < nr_of_triangles; ++tria_index)
        {
        if (!marked[tria_index])
          {
          ++current_shell;
          std::queue<uint32_t> qu;
          qu.push(tria_index);
          marked[tria_index] = true;
          while (!qu.empty())
            {
            uint32_t current_tria = qu.front();
            qu.pop();
            sh[current_tria] = current_shell - 1;
            auto neighbs = tr_neigh[current_tria];
            for (uint32_t n : neighbs)
              {
              if (!marked[n])
                {
                if (manifold)
                  {
                  uint32_t v0, v1;
                  edge_between_triangles(v0, v1, triangles[current_tria], triangles[n]);
                  if (is_manifold_edge(v0, v1, adj_list))
                    {
                    qu.push(n);
                    marked[n] = true;
                    }
                  }
                else
                  {
                  qu.push(n);
                  marked[n] = true;
                  }
                }
              }
            }
          }
        }
      return sh;
      }
    }


  inline std::vector<uint32_t> shells(const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, shell_connectivity conn)
    {
    using namespace details;
    switch (conn)
      {
      case shell_connectivity::edge: return _shells(nr_of_vertices, triangles, nr_of_triangles, false);
      case shell_connectivity::manifold: return _shells(nr_of_vertices, triangles, nr_of_triangles, true);
      case shell_connectivity::vertex:
      {
      std::vector<uint32_t> sh(nr_of_triangles, 0);
      auto adj_list = adjacency_list(nr_of_vertices, triangles, nr_of_triangles);
      std::vector<bool> marked(sh.size(), false);
      uint32_t current_shell = 0;
      for (uint32_t tria_index = 0; tria_index < nr_of_triangles; ++tria_index)
        {
        if (!marked[tria_index])
          {
          ++current_shell;
          std::queue<uint32_t> qu;
          qu.push(tria_index);
          marked[tria_index] = true;
          while (!qu.empty())
            {
            uint32_t current_tria = qu.front();
            auto current_tria_idces = triangles[current_tria];
            qu.pop();
            sh[current_tria] = current_shell - 1;
            std::set<uint32_t> neighbs;
            neighbs.insert(adj_list.begin(current_tria_idces[0]), adj_list.end(current_tria_idces[0]));
            for (uint32_t idx = 1; idx < 3; ++idx)
              {
              neighbs.insert(adj_list.begin(current_tria_idces[idx]), adj_list.end(current_tria_idces[idx]));
              }
            for (uint32_t n : neighbs)
              {
              if (!marked[n])
                {
                qu.push(n);
                marked[n] = true;
                }
              }
            }
          }
        }
      return sh;
      }
      default: return _shells(nr_of_vertices, triangles, nr_of_triangles, true);
      }
    }

  inline trivial_ear::trivial_ear() : vertices(nullptr), adj_list(nullptr), v0(-1), v1(-1), v2(-1) {}

  inline trivial_ear::trivial_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list) :
    v0(p0), v1(p1), v2(p2), vertices(p_vertices), adj_list(p_adj_list)
    {
    n = normalize(cross(vertices[v2] - vertices[v1], vertices[v0] - vertices[v1]));
    compute_quality();
    compute_angle(p_triangles);
    }

  inline void trivial_ear::compute_quality()
    {
    auto p0 = vertices[v0];
    auto p1 = vertices[v1];
    auto p2 = vertices[v2];
    float e1 = distance_sqr(p0, p1);
    float e2 = distance_sqr(p0, p2);
    float e3 = distance_sqr(p1, p2);
    float area = length(cross(p2 - p1, p0 - p1)) * 0.5f;
    quality = area / (e1 + e2 + e3);
    }

  inline void trivial_ear::compute_angle(const vec3<uint32_t>* triangles)
    {
    auto point1 = vertices[v2] - vertices[v1];
    auto point2 = vertices[v0] - vertices[v1];
    auto tmp1 = point1 * length(point2);
    auto tmp2 = point2 * length(point1);
    angle_rad = 2.f * std::atan(length(tmp1 - tmp2) / length(tmp1 + tmp2));
    auto tria_idcs = triangle_indices_from_edge(v0, v1, *adj_list);
    if (tria_idcs.empty())
      return;
    uint32_t t = tria_idcs.front();
    auto tria = triangles[t];
    auto tria_n = normalize(cross(vertices[tria[2]] - vertices[tria[1]], vertices[tria[0]] - vertices[tria[1]]));
    float flip_angle = dot(n, tria_n);
    if (flip_angle < 0.f)
      angle_rad = 2.f*3.141592653589793238462643383f - angle_rad;
    }

  inline bool trivial_ear::is_concave() const
    {
    return angle_rad > 3.141592653589793238462643383f;
    }

  inline bool trivial_ear::operator < (const trivial_ear& other) const
    {
    if (is_concave() && !other.is_concave())
      return true;
    if (!is_concave() && other.is_concave())
      return false;
    return quality < other.quality;
    }

  inline float dihedral_angle(const vec3<float>& u, const vec3<float>& v, const vec3<float>& a, const vec3<float>& b)
    {
    auto n0 = normalize(cross(v - u, a - u));
    auto n1 = normalize(cross(u - v, b - v));
    float dt = dot(n0, n1);
    if (dt > 1.f)
      return 0.f;
    else if (dt < -1.f)
      return 3.141592653589793238462643383f;
    return acos(dt);
    }

  inline minimum_weight_ear::minimum_weight_ear() : trivial_ear()
    {}

  inline minimum_weight_ear::minimum_weight_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list) :
    trivial_ear(p0, p1, p2, p_vertices, p_triangles, p_adj_list)
    {
    compute_quality(p_triangles);
    }

  inline void minimum_weight_ear::compute_quality(const vec3<uint32_t>* triangles)
    {
    aspect_ratio = quality;
    dihedral_rad = 0.f;
    auto tria_idcs_a = triangle_indices_from_edge(v0, v1, *adj_list);
    if (tria_idcs_a.empty())
      return;
    auto tria_idcs_b = triangle_indices_from_edge(v1, v2, *adj_list);
    if (tria_idcs_b.empty())
      return;
    uint32_t t = tria_idcs_a.front();
    auto tria_a = triangles[t];
    size_t a = tria_a[0];
    if (a == v0 || a == v1)
      {
      a = tria_a[1];
      if (a == v0 || a == v1)
        a = tria_a[2];
      }
    t = tria_idcs_b.front();
    auto tria_b = triangles[t];
    size_t b = tria_b[0];
    if (b == v2 || b == v1)
      {
      b = tria_b[1];
      if (b == v2 || b == v1)
        b = tria_b[2];
      }
    double dihedral_a = dihedral_angle(vertices[v0], vertices[v1], vertices[v2], vertices[a]);
    double dihedral_b = dihedral_angle(vertices[v1], vertices[v2], vertices[v0], vertices[b]);
    dihedral_rad = std::max(dihedral_a, dihedral_b);
    }

  inline bool minimum_weight_ear::operator < (const minimum_weight_ear& other) const
    {
    if (is_concave() && !other.is_concave())
      return true;
    if (!is_concave() && other.is_concave())
      return false;
    return (aspect_ratio - (dihedral_rad / 3.141592653589793238462643383f)*0.1f) < (other.aspect_ratio - (other.dihedral_rad / 3.141592653589793238462643383f)*0.1f);
    }

  template <class Ear>
  void fill_hole_ear(std::vector<vec3<uint32_t>>& triangles, const vec3<float>* vertices, mutable_adjacency_list& mal, const std::vector<uint32_t>& hole)
    {
    size_t sz = hole.size();
    if (sz < 3)
      return;

    std::vector<Ear> ear_heap;
    ear_heap.reserve(sz);
    for (size_t v1 = 0; v1 < sz; ++v1)
      {
      size_t v0 = (v1 + (sz - 1)) % sz;
      size_t v2 = (v1 + 1) % sz;
      Ear new_ear(hole[v0], hole[v1], hole[v2], vertices, triangles.data(), &mal);
      ear_heap.push_back(new_ear);
      }
    std::make_heap(ear_heap.begin(), ear_heap.end());
    size_t cnt = sz;
    while (cnt > 2 && !ear_heap.empty())
      {
      Ear best_ear = ear_heap.front();
      std::pop_heap(ear_heap.begin(), ear_heap.end());
      ear_heap.pop_back();

      auto trias02 = triangle_indices_from_edge(best_ear.v0, best_ear.v2, mal);
      auto trias01 = triangle_indices_from_edge(best_ear.v0, best_ear.v1, mal);
      auto trias12 = triangle_indices_from_edge(best_ear.v1, best_ear.v2, mal);
      if (trias02.size() < 2 && trias01.size() == 1 && trias12.size() == 1)
        {
        uint32_t id = (uint32_t)triangles.size();
        triangles.emplace_back(best_ear.v0, best_ear.v1, best_ear.v2);
        mal.add_triangle_to_vertex(best_ear.v0, id);
        mal.add_triangle_to_vertex(best_ear.v1, id);
        mal.add_triangle_to_vertex(best_ear.v2, id);
        --cnt;
        auto neighbour_vertices_v0 = one_ring_vertices_from_vertex(best_ear.v0, mal, triangles.data());
        for (auto nv0 : neighbour_vertices_v0)
          {
          if (nv0 == best_ear.v2)
            continue;
          if (triangle_indices_from_edge(best_ear.v0, nv0, mal).size() == 1)
            {
            Ear new_ear(nv0, best_ear.v0, best_ear.v2, vertices, triangles.data(), &mal);
            ear_heap.push_back(new_ear);
            std::push_heap(ear_heap.begin(), ear_heap.end());
            }
          }
        auto neighbour_vertices_v2 = one_ring_vertices_from_vertex(best_ear.v2, mal, triangles.data());
        for (auto nv2 : neighbour_vertices_v2)
          {
          if (nv2 == best_ear.v0)
            continue;
          if (triangle_indices_from_edge(best_ear.v2, nv2, mal).size() == 1)
            {
            Ear new_ear(best_ear.v0, best_ear.v2, nv2, vertices, triangles.data(), &mal);
            ear_heap.push_back(new_ear);
            std::push_heap(ear_heap.begin(), ear_heap.end());
            }
          }
        }
      }
    }

  template <class Ear>
  void fill_holes(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, uint32_t max_holes_size)
    {
    mutable_adjacency_list adj_list((uint32_t)vertices.size(), triangles.data(), (uint32_t)triangles.size());
    std::vector<bool> vertex_treated(vertices.size(), false);
    std::vector<std::vector<uint32_t>> holes;
    for (uint32_t v = 0; v < (uint32_t)vertex_treated.size(); ++v)
      {
      if (!vertex_treated[v])
        {
        vertex_treated[v] = true;
        if (is_boundary_vertex(v, adj_list, triangles.data()))
          {
          std::vector<uint32_t> hole;
          hole.push_back(v);
          std::queue<uint32_t> qu;
          auto neighbouring_vertices = one_ring_vertices_from_vertex(v, adj_list, triangles.data());
          for (auto v2 : neighbouring_vertices)
            {
            if (!vertex_treated[v2] && is_boundary_edge(v, v2, adj_list))
              {
              qu.push(v2);
              break;
              }
            }
          while (!qu.empty())
            {
            uint32_t current_vertex = qu.front();
            assert(!vertex_treated[current_vertex]);
            hole.push_back(current_vertex);
            vertex_treated[current_vertex] = true;
            qu.pop();
            neighbouring_vertices = one_ring_vertices_from_vertex(current_vertex, adj_list, triangles.data());
            for (auto v2 : neighbouring_vertices)
              {
              if (!vertex_treated[v2] && is_boundary_edge(current_vertex, v2, adj_list))
                {
                qu.push(v2);
                break;
                }
              }
            }
          bool valid_hole = hole.size() > 2;
          if (valid_hole)
            {
            neighbouring_vertices = one_ring_vertices_from_vertex(v, adj_list, triangles.data());
            bool found_last = false;
            for (auto v2 : neighbouring_vertices)
              {
              if (v2 == hole.back())
                {
                found_last = true;
                break;
                }
              }
            valid_hole = found_last;
            }
          if (valid_hole)
            holes.push_back(hole);
          }
        }
      }

    for (auto& hole : holes)
      {
      if (hole.size() < max_holes_size)
        {
        uint32_t v0 = hole[0];
        uint32_t v1 = hole[1];
        auto tri = triangle_indices_from_edge(v0, v1, adj_list);
        bool flip = false;
        auto tria = triangles[tri[0]]; // since these are boundary edges there will be exactly one triangle
        for (int j = 0; j < 3; ++j)
          {
          if (tria[j] == v0 && tria[(j + 1) % 3] == v1)
            flip = true;
          }
        if (flip)
          std::reverse(hole.begin(), hole.end());
        fill_hole_ear<Ear>(triangles, vertices.data(), adj_list, hole);
        }
      }
    }

  namespace details
    {

    inline void _CreateFilter(std::vector<std::vector<uint32_t>>& io_filter, const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles)
      {
      std::vector<std::vector<uint32_t>> filter(vertices.size());
      for (size_t i = 0; i < triangles.size(); ++i)
        {
        uint32_t ind1 = triangles[i][0];
        uint32_t ind2 = triangles[i][1];
        uint32_t ind3 = triangles[i][2];
        filter[ind1].push_back(ind2);
        filter[ind2].push_back(ind3);
        filter[ind3].push_back(ind1);
        }
      io_filter.swap(filter);
      }

    inline void _ApplyFilter(std::vector<vec3<float>>& o_vertices, std::vector<vec3<float>>& i_vertices, const std::vector<std::vector<uint32_t>>& i_filter, float filter_value)
      {
      if (!filter_value)
        {
        o_vertices.swap(i_vertices);
        return;
        }
      parallel_for(size_t(0), i_vertices.size(), [&](size_t j)
        {
        size_t valence = i_filter[j].size();
        for (size_t k = 0; k < 3; ++k)
          {
          o_vertices[j][k] = i_vertices[j][k] * (1.f - filter_value);
          for (size_t v = 0; v < valence; ++v)
            o_vertices[j][k] += filter_value / float(valence)*i_vertices[i_filter[j][v]][k];
          }
        });
      }

    inline void _ApplyLocalFilter(std::vector<vec3<float>>& o_vertices, std::vector<vec3<float>>& i_vertices, const std::vector<uint32_t>& vertex_indices, const std::vector<std::vector<uint32_t>>& i_filter, float filter_value)
      {
      if (!filter_value)
        {
        o_vertices.swap(i_vertices);
        return;
        }
      parallel_for(size_t(0), vertex_indices.size(), [&](size_t j)
        {
        size_t vertex_index = vertex_indices[j];
        size_t valence = i_filter[vertex_index].size();
        if (valence)
          for (size_t k = 0; k < 3; ++k)
            {
            o_vertices[vertex_index][k] = i_vertices[vertex_index][k] * (1.f - filter_value);
            for (size_t v = 0; v < valence; ++v)
              o_vertices[vertex_index][k] += filter_value / float(valence)*i_vertices[i_filter[vertex_index][v]][k];
            }
        });
      }
    }

  inline void smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, uint32_t iterations, float lambda, float mu)
    {
    using namespace details;
    std::vector<std::vector<uint32_t>> filter;
    _CreateFilter(filter, vertices, triangles);
    std::vector<vec3<float>> new_vertices(vertices);
    for (size_t i = 0; i < iterations; ++i)
      {
      _ApplyFilter(new_vertices, vertices, filter, lambda);
      _ApplyFilter(vertices, new_vertices, filter, mu);
      }
    }

  inline void local_smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& vertex_indices, uint32_t iterations, float lambda, float mu)
    {
    using namespace details;
    std::vector<std::vector<uint32_t>> filter;
    _CreateFilter(filter, vertices, triangles);
    std::vector<vec3<float>> new_vertices(vertices);
    for (size_t i = 0; i < iterations; ++i)
      {
      _ApplyLocalFilter(new_vertices, vertices, vertex_indices, filter, lambda);
      _ApplyLocalFilter(vertices, new_vertices, vertex_indices, filter, mu);
      }
    }

  inline void dyadic_subdivide(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles)
    {
    const uint32_t nr_of_vertices = vertices.size();
    const uint32_t nr_of_triangles = triangles.size();
    vertices.reserve(nr_of_vertices * 2 + nr_of_triangles + 10);

    std::unordered_map<uint64_t, uint32_t> treated(nr_of_vertices);

    std::sort(triangles.begin(), triangles.end()); // necessary for undo redo after scambling triangles because of bvh

    for (const auto& tria : triangles)
      {
      for (int j = 0; j < 3; ++j)
        {
        uint32_t v0 = tria[j];
        uint32_t v1 = tria[(j + 1) % 3];
        if (v0 > v1)
          std::swap(v0, v1);
        uint64_t key = (((uint64_t)v0) << 32) | (uint64_t)v1;
        if (treated.find(key) == treated.end())
          {
          treated.insert(std::make_pair(key, (uint32_t)vertices.size()));
          vertices.push_back((vertices[v0] + vertices[v1]) / 2.f);
          }
        }
      }

    triangles.reserve(nr_of_triangles * 4);
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      auto& tria = triangles[t];
      uint32_t edge_points[3];
      for (int j = 0; j < 3; ++j)
        {
        uint32_t v0 = tria[j];
        uint32_t v1 = tria[(j + 1) % 3];
        if (v0 > v1)
          std::swap(v0, v1);
        uint64_t key = (((uint64_t)v0) << 32) | (uint64_t)v1;
        edge_points[j] = treated.find(key)->second;
        }
      triangles.push_back(vec3<uint32_t>(tria[0], edge_points[0], edge_points[2]));
      triangles.push_back(vec3<uint32_t>(edge_points[0], tria[1], edge_points[1]));
      triangles.push_back(vec3<uint32_t>(edge_points[2], edge_points[1], tria[2]));
      tria[0] = edge_points[0];
      tria[1] = edge_points[1];
      tria[2] = edge_points[2];
      }
    }

  } // namespace jtk