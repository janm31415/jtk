/*
   Do this:
      #define JTK_GEOMETRY_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_GEOMETRY_IMPLEMENTATION
   #include "jtk/geometry.h"
 */

#ifndef JTK_GEOMETRY_H
#define JTK_GEOMETRY_H

#include "vec.h"
#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#ifndef JTKGDEF
#ifdef JTK_GEOMETRY_STATIC
#define JTKGDEF static
#else
#define JTKGDEF extern
#endif
#endif

#ifndef JTKGINLINE
#define JTKGINLINE inline 
#endif

namespace jtk
  {

  /////////////////////////////////////////////////////////////////////////
  // interfaces
  /////////////////////////////////////////////////////////////////////////

  JTKGDEF bool read_stl(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  JTKGDEF bool read_stl_ascii(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  JTKGDEF bool write_stl(const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const vec3<float>* triangle_normals, const unsigned short* attributes, const char* filename);
  JTKGDEF bool read_obj(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  JTKGDEF bool read_obj(std::string& mtl_filename, std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv, const char* filename);
  JTKGDEF bool read_texture_filename_from_mtl(std::string& texture_file, const char* filename);
  JTKGDEF bool read_off(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename);
  JTKGDEF bool write_off(uint32_t nr_of_vertices, const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const char* filename);

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts);
  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<float>>& normals);
  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<uint32_t>& clrs);
  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs);
  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<uint32_t>>& triangles);
  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<uint32_t>& pts_colors, const std::vector<vec3<uint32_t>>& triangles);
  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<jtk::vec2<float>>>& uv);

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

  JTKGDEF void delete_triangles(std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& triangles_to_delete);

  JTKGDEF void remove_free_vertices(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles);

  JTKGDEF void remove_free_vertices(std::vector<vec3<float>>& vertices, std::vector<uint32_t>& vertex_colors, std::vector<vec3<uint32_t>>& triangles);

  JTKGDEF bool edge_swap(uint32_t v, uint32_t v2, std::vector<jtk::vec3<uint32_t>>& triangles, jtk::mutable_adjacency_list& adj_list);

  JTKGDEF std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  JTKGDEF std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  JTKGDEF std::vector<std::vector<uint32_t>> ordered_one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles, bool oriented = true);

  JTKGDEF void compute_triangle_normals(std::vector<vec3<float>>& triangle_normals, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles);

  JTKGDEF void compute_vertex_normals(std::vector<vec3<float>>& vertex_normals, const vec3<float>* triangle_normals, const vec3<float>* vertices, const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles);

  JTKGDEF vec3<float> normal(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles);

  JTKGDEF vec3<float> normal(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  JTKGDEF float area(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles);

  JTKGDEF std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list);

  JTKGDEF std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list);

  JTKGDEF float signed_volume(const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles);

  JTKGDEF bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  JTKGDEF bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* triangle_normals, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  JTKGDEF bool vertex_on_sharp_edge(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle = 0.86602540f);

  JTKGDEF std::vector<std::vector<uint32_t>> triangle_neighbour_indices(const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, const adjacency_list& adj_list);

  // returns false if tria1 and tria2 share three or no edges (in the latter case e contains std::numeric_limits<size_t>::max() twice)
  // otherwise returns the shared edge oriented such as it appears in tria1
  JTKGDEF bool edge_between_triangles(uint32_t& v0, uint32_t& v1, const vec3<uint32_t>& tria1, const vec3<uint32_t>& tria2);

  JTKGDEF bool is_manifold_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list);

  JTKGDEF bool is_boundary_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list);

  JTKGDEF bool is_boundary_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list);

  JTKGDEF bool is_boundary_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  JTKGDEF bool is_boundary_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  enum class shell_connectivity { vertex, edge, manifold };
  JTKGDEF std::vector<uint32_t> shells(const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, shell_connectivity conn = shell_connectivity::edge);

  JTKGDEF std::vector<std::vector<uint32_t>> holes(const adjacency_list& adj_list, const vec3<uint32_t>* triangles);

  JTKGDEF std::vector<std::vector<uint32_t>> holes(const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles);

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

  JTKGDEF float dihedral_angle(const vec3<float>& u, const vec3<float>& v, const vec3<float>& a, const vec3<float>& b);

  struct minimum_weight_ear : public trivial_ear
    {
    minimum_weight_ear();
    minimum_weight_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list);
    virtual void compute_quality(const vec3<uint32_t>* triangles);
    virtual bool operator < (const minimum_weight_ear& other) const;

    float aspect_ratio;
    float dihedral_rad;
    };

  template <class Ear>
  void fill_hole_ear(std::vector<vec3<uint32_t>>& triangles, const vec3<float>* vertices, mutable_adjacency_list& mal, const std::vector<uint32_t>& hole);
  template <class Ear>
  void fill_holes(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, uint32_t max_holes_size);

  JTKGDEF void smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, uint32_t iterations = 100, float lambda = 0.33f, float mu = -0.34f);
  JTKGDEF void local_smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& vertex_indices, uint32_t iterations = 100, float lambda = 0.33f, float mu = -0.34f);
  JTKGDEF void dyadic_subdivide(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles);
  JTKGDEF void dyadic_subdivide_uv_map(std::vector<jtk::vec3<jtk::vec2<float>>>& triangle_uv);
  JTKGDEF void dyadic_subdivide_triangle_indices_vector(std::vector<uint32_t>& triangle_indices, uint32_t original_number_of_triangles);
  JTKGDEF void undo_dyadic_subdivide(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles);
  JTKGDEF void undo_dyadic_subdivide_uv_map(std::vector<jtk::vec3<jtk::vec2<float>>>& triangle_uv);
  JTKGDEF void butterfly(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles);

  template <class TDistanceFunction, class TValidValue>
  void marching_cubes(
    std::vector<vec3<float>>& vertices,
    std::vector<vec3<uint32_t>>& triangles,
    const boundingbox3d<float>& bb,
    const uint32_t width, const uint32_t height, const uint32_t depth,
    float isovalue,
    TDistanceFunction fun,
    TValidValue valid_value);


  JTKGDEF bool stitch_points(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, float distance_tolerance);



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

    auto rm_iter = indices_to_delete.begin();
    auto end = indices_to_delete.end();
    std::size_t current_index = 0;

    const auto pred = [&](const TType&) {
      // any more to remove?
      if (rm_iter == end) { return false; }
      // is this one specified?
      if (*rm_iter == current_index++) { return ++rm_iter, true; }
      return false;
      };

    vec.erase(std::remove_if(vec.begin(), vec.end(), pred), vec.end());
    /*
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
    */
    }

  namespace marching_cubes_details
    {
    // For any edge, if one vertex is inside of the surface and the other is outside of the surface
    //  then the edge intersects the surface
    // For each of the 8 vertices of the cube can be two possible states : either inside or outside of the surface
    // For any cube the are 2^8=256 possible sets of vertex states
    // This table lists the edges intersected by the surface for all 256 possible vertex states
    // There are 12 edges.  For each entry in the table, if edge #n is intersected, then bit #n is set to 1

    static int aiCubeEdgeFlags[256] =
      {
      0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
      0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
      0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
      0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
      0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
      0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
      0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
      0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
      0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
      0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
      0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
      0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
      0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
      0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
      0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
      0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
      };

    //  For each of the possible vertex states listed in aiCubeEdgeFlags there is a specific triangulation
    //  of the edge intersection points.  a2iTriangleConnectionTable lists all of them in the form of
    //  0-5 edge triples with the list terminated by the invalid value -1.
    //  For example: a2iTriangleConnectionTable[3] list the 2 triangles formed when corner[0]
    //  and corner[1] are inside of the surface, but the rest of the cube is not.
    //
    //  I found this table in an example program someone wrote long ago.  It was probably generated by hand

    static int a2iTriangleConnectionTable[256][16] =
      {
          { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
          { 8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1 },
          { 3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1 },
          { 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
          { 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1 },
          { 9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1 },
          { 10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1 },
          { 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
          { 5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1 },
          { 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1 },
          { 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
          { 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
          { 11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1 },
          { 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1 },
          { 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1 },
          { 11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1 },
          { 2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1 },
          { 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
          { 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
          { 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
          { 6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1 },
          { 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1 },
          { 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
          { 3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
          { 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1 },
          { 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1 },
          { 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
          { 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1 },
          { 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1 },
          { 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
          { 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
          { 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1 },
          { 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1 },
          { 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
          { 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1 },
          { 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1 },
          { 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1 },
          { 3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
          { 6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1 },
          { 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1 },
          { 10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
          { 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1 },
          { 7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
          { 7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
          { 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
          { 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1 },
          { 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1 },
          { 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1 },
          { 0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1 },
          { 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1 },
          { 7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1 },
          { 10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1 },
          { 0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1 },
          { 7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
          { 6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1 },
          { 6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1 },
          { 4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1 },
          { 10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1 },
          { 8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1 },
          { 1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1 },
          { 10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1 },
          { 10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1 },
          { 9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1 },
          { 7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1 },
          { 3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1 },
          { 7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1 },
          { 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1 },
          { 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1 },
          { 9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1 },
          { 1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1 },
          { 4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1 },
          { 7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1 },
          { 6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1 },
          { 0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1 },
          { 6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1 },
          { 0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1 },
          { 11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1 },
          { 6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1 },
          { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1 },
          { 9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1 },
          { 1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1 },
          { 10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1 },
          { 0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1 },
          { 10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1 },
          { 11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1 },
          { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1 },
          { 7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1 },
          { 2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1 },
          { 9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1 },
          { 9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1 },
          { 1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
          { 5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1 },
          { 0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1 },
          { 10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1 },
          { 2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1 },
          { 0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1 },
          { 0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1 },
          { 9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1 },
          { 5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1 },
          { 5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1 },
          { 8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1 },
          { 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1 },
          { 1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1 },
          { 3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1 },
          { 4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1 },
          { 9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1 },
          { 11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
          { 11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1 },
          { 2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1 },
          { 9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1 },
          { 3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1 },
          { 1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1 },
          { 4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1 },
          { 0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
          { 9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1 },
          { 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
          { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }
      };

    //a2fVertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
    static const float a2fVertexOffset[8][3] =
      {
          { 0.f, 0.f, 0.f }, { 1.f, 0.f, 0.f }, { 1.f, 1.f, 0.f }, { 0.f, 1.f, 0.f },
          { 0.f, 0.f, 1.f }, { 1.f, 0.f, 1.f }, { 1.f, 1.f, 1.f }, { 0.f, 1.f, 1.f }
      };


    static const int a2fVoxelOffset[8][3] =
      {
          { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
          { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 }
      };


    //a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
    static const int a2iEdgeConnection[12][2] =
      {
          { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
          { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
          { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
      };

    //a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
    static const float a2fEdgeDirection[12][3] =
      {
          { 1.f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, { -1.f, 0.f, 0.f }, { 0.f, -1.f, 0.f },
          { 1.f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, { -1.f, 0.f, 0.f }, { 0.f, -1.f, 0.f },
          { 0.f, 0.f, 1.f }, { 0.f, 0.f, 1.f }, { 0.f, 0.f, 1.f }, { 0.f, 0.f, 1.f }
      };

    //fGetOffset finds the approximate point of intersection of the surface
    // between two points with the values fValue1 and fValue2
    inline float fGetOffset(float fValue1, float fValue2, float fValueDesired)
      {
      float fDelta = fValue2 - fValue1;

      if (fDelta == 0.f)
        {
        return 0.5f;
        }
      return (fValueDesired - fValue1) / fDelta;
      }

    struct vec3floathash
      {
      size_t operator()(const vec3<float>& point) const
        {
        std::hash<float> hash_fn;
        return (hash_fn(point[0]) * 73856093) ^ (hash_fn(point[1]) * 19349663) ^ (hash_fn(point[2]) * 83492791);
        }
      };
    }

  template <class TDistanceFunction, class TValidValue>
  void marching_cubes(
    std::vector<vec3<float>>& vertices,
    std::vector<vec3<uint32_t>>& triangles,
    const boundingbox3d<float>& bb,
    const uint32_t width, const uint32_t height, const uint32_t depth,
    float isovalue,
    TDistanceFunction fun,
    TValidValue valid_value)
    {
    assert(width);
    assert(height);
    assert(depth);
    assert(bb.max[0] > bb.min[0]);
    assert(bb.max[1] > bb.min[1]);
    assert(bb.max[2] > bb.min[2]);

    using namespace marching_cubes_details;

    auto less_fie = [](const vec3<float>& left, const vec3<float>& right)
      {
      return left[0] == right[0] ? (left[1] == right[1] ? left[2] < right[2] : left[1] < right[1]) : left[0] < right[0];
      };

    const float delta_x = (bb.max[0] - bb.min[0]) / (float)width;
    const float delta_y = (bb.max[1] - bb.min[1]) / (float)height;
    const float delta_z = (bb.max[2] - bb.min[2]) / (float)depth;

    vertices.clear();
    triangles.clear();

    std::vector<std::vector<vec3<float>>> pts(depth - 1);

    parallel_for(int(0), int(depth - 1), [&](int z)
      //for (int z = 0; z < int(depth - 1); ++z)
      {
      for (int y = 0; y<int(height - 1); ++y)
        {
        for (int x = 0; x<int(width - 1); ++x)
          {
          int iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
          float fOffset;
          float afCubeValue[8];
          vec3<float> asEdgeVertex[12];
          for (int v = 0; v < 8; ++v)
            {
            int dx = a2fVoxelOffset[v][0];
            int dy = a2fVoxelOffset[v][1];
            int dz = a2fVoxelOffset[v][2];
            int new_x = x + dx;
            int new_y = y + dy;
            int new_z = z + dz;

            afCubeValue[v] = (float)fun(new_x * delta_x + bb.min[0], new_y * delta_x + bb.min[1], new_z * delta_z + bb.min[2]);
            }
          //Find which vertices are inside of the surface and which are outside
          iFlagIndex = 0;
          for (iVertexTest = 0; iVertexTest < 8; ++iVertexTest)
            {
            if (!valid_value(afCubeValue[iVertexTest]))
              {
              iFlagIndex = 0;
              break;
              }
            else if ((float)afCubeValue[iVertexTest] <= isovalue)
              iFlagIndex |= 1 << iVertexTest;
            }
          //Find which edges are intersected by the surface
          iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

          //If the cube is entirely inside or outside of the surface, then there will be no intersections
          if (iEdgeFlags == 0)
            {
            continue;
            }

          //Find the point of intersection of the surface with each edge
          //Then find the normal to the surface at those points
          for (iEdge = 0; iEdge < 12; ++iEdge)
            {
            //if there is an intersection on this edge
            if (iEdgeFlags & (1 << iEdge))
              {
              vec3<float> corner_pt((float)x, (float)y, (float)z);
              vec3<float> pt0((float)x + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][0]), (float)y + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][1]), (float)z + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][2]));
              vec3<float> pt1((float)x + (a2fVertexOffset[a2iEdgeConnection[iEdge][1]][0]), (float)y + (a2fVertexOffset[a2iEdgeConnection[iEdge][1]][1]), (float)z + (a2fVertexOffset[a2iEdgeConnection[iEdge][1]][2]));
              float val0 = (float)afCubeValue[a2iEdgeConnection[iEdge][0]];
              float val1 = (float)afCubeValue[a2iEdgeConnection[iEdge][1]];
              if (less_fie(pt1, pt0))
                {
                vec3<float> tmp(pt1);
                pt1 = pt0;
                pt0 = tmp;
                std::swap(val0, val1);
                }
              fOffset = fGetOffset(val0, val1, isovalue);
              vec3<float> p = pt0 + (pt1 - pt0) * fOffset;
              asEdgeVertex[iEdge][0] = p[0] * delta_x + bb.min[0];
              asEdgeVertex[iEdge][1] = p[1] * delta_y + bb.min[1];
              asEdgeVertex[iEdge][2] = p[2] * delta_z + bb.min[2];
              }
            }

          //Draw the triangles that were found.  There can be up to five per cube
          for (iTriangle = 0; iTriangle < 5; ++iTriangle)
            {
            if (a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle] < 0)
              break;

            vec3<float> vert[3];
            for (iCorner = 0; iCorner < 3; ++iCorner)
              {
              iVertex = a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle + iCorner];
              vert[iCorner] = asEdgeVertex[iVertex];
              }

            pts[z].push_back(vert[0]);
            pts[z].push_back(vert[2]);
            pts[z].push_back(vert[1]);
            }
          }
        }
      });

    std::unordered_map<vec3<float>, uint32_t, vec3floathash> point_map;

    for (auto it = pts.begin(); it != pts.end(); ++it)
      {
      for (size_t i = 0; i < it->size() / 3; ++i)
        {
        vec3<float> v0 = (*it)[i * 3];
        vec3<float> v1 = (*it)[i * 3 + 1];
        vec3<float> v2 = (*it)[i * 3 + 2];

        uint32_t V0, V1, V2;
        auto it = point_map.find(v0);
        if (it == point_map.end())
          {
          V0 = (uint32_t)vertices.size();
          vertices.push_back(v0);
          point_map[v0] = V0;
          }
        else
          V0 = it->second;

        it = point_map.find(v1);
        if (it == point_map.end())
          {
          V1 = (uint32_t)vertices.size();
          vertices.push_back(v1);
          point_map[v1] = V1;
          }
        else
          V1 = it->second;

        it = point_map.find(v2);
        if (it == point_map.end())
          {
          V2 = (uint32_t)vertices.size();
          vertices.push_back(v2);
          point_map[v2] = V2;
          }
        else
          V2 = it->second;

        triangles.emplace_back(V0, V1, V2);
        }
      }
    }


  JTKGINLINE adjacency_list::adjacency_list() : _nr_of_vertices(0)
    {

    }

  JTKGINLINE adjacency_list::adjacency_list(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles) : _nr_of_vertices(nr_of_vertices)
    {
    build(nr_of_vertices, triangles, nr_of_triangles);
    }

  JTKGINLINE adjacency_list::~adjacency_list()
    {

    }

  JTKGINLINE void adjacency_list::build(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
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

  JTKGINLINE adjacency_list::const_iterator adjacency_list::begin(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list.begin() + offset_per_vertex[vertex_index];
    }

  JTKGINLINE adjacency_list::const_iterator adjacency_list::end(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list.begin() + offset_per_vertex[vertex_index + 1];
    }

  JTKGINLINE uint32_t adjacency_list::size() const
    {
    return _nr_of_vertices;
    }

  JTKGINLINE bool adjacency_list::empty() const
    {
    return _nr_of_vertices == 0;
    }

  JTKGINLINE mutable_adjacency_list::mutable_adjacency_list()
    {

    }

  JTKGINLINE mutable_adjacency_list::mutable_adjacency_list(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
    {
    build(nr_of_vertices, triangles, nr_of_triangles);
    }

  JTKGINLINE mutable_adjacency_list::~mutable_adjacency_list()
    {

    }

  JTKGINLINE void mutable_adjacency_list::build(uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, uint32_t nr_of_triangles)
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

  JTKGINLINE void mutable_adjacency_list::add_triangle_to_vertex(uint32_t vertex_id, uint32_t triangle_id)
    {
    auto first = triangles_per_vertex_list[vertex_id].begin();
    auto last = triangles_per_vertex_list[vertex_id].end();
    auto it = std::lower_bound(first, last, triangle_id);
    triangles_per_vertex_list[vertex_id].insert(it, triangle_id);
    }

  JTKGINLINE void mutable_adjacency_list::remove_triangle_from_vertex(uint32_t vertex_id, uint32_t triangle_id)
    {
    auto first = triangles_per_vertex_list[vertex_id].begin();
    auto last = triangles_per_vertex_list[vertex_id].end();
    auto it = std::lower_bound(first, last, triangle_id);
    if (it != last && *it == triangle_id)
      {
      triangles_per_vertex_list[vertex_id].erase(it);
      }
    }

  JTKGINLINE void mutable_adjacency_list::add_triangles_to_vertex(uint32_t vertex_index, const std::vector<uint32_t>& triangle_indices)
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

  JTKGINLINE void mutable_adjacency_list::remove_triangles_from_vertex(uint32_t vertex_index, const std::vector<uint32_t>& triangle_indices)
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

  JTKGINLINE mutable_adjacency_list::const_iterator mutable_adjacency_list::begin(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list[vertex_index].begin();
    }

  JTKGINLINE mutable_adjacency_list::const_iterator mutable_adjacency_list::end(uint32_t vertex_index) const
    {
    return triangles_per_vertex_list[vertex_index].end();
    }

  JTKGINLINE uint32_t mutable_adjacency_list::size(uint32_t vertex_index) const
    {
    return (uint32_t)triangles_per_vertex_list[vertex_index].size();
    }

  JTKGINLINE uint32_t mutable_adjacency_list::size() const
    {
    return (uint32_t)triangles_per_vertex_list.size();
    }

  JTKGINLINE bool mutable_adjacency_list::empty() const
    {
    return triangles_per_vertex_list.empty();
    }

  JTKGINLINE trivial_ear::trivial_ear() : vertices(nullptr), adj_list(nullptr), v0((uint32_t)-1), v1((uint32_t)-1), v2((uint32_t)-1) {}

  JTKGINLINE trivial_ear::trivial_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list) :
    v0(p0), v1(p1), v2(p2), vertices(p_vertices), adj_list(p_adj_list)
    {
    n = normalize(cross(vertices[v2] - vertices[v1], vertices[v0] - vertices[v1]));
    compute_quality();
    compute_angle(p_triangles);
    }

  JTKGINLINE void trivial_ear::compute_quality()
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

  JTKGINLINE void trivial_ear::compute_angle(const vec3<uint32_t>* triangles)
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
      angle_rad = 2.f * 3.141592653589793238462643383f - angle_rad;
    }

  JTKGINLINE bool trivial_ear::is_concave() const
    {
    return angle_rad > 3.141592653589793238462643383f;
    }

  JTKGINLINE bool trivial_ear::operator < (const trivial_ear& other) const
    {
    if (is_concave() && !other.is_concave())
      return true;
    if (!is_concave() && other.is_concave())
      return false;
    return quality < other.quality;
    }

  JTKGINLINE minimum_weight_ear::minimum_weight_ear() : trivial_ear()
    {}

  JTKGINLINE minimum_weight_ear::minimum_weight_ear(uint32_t p0, uint32_t p1, uint32_t p2, const vec3<float>* p_vertices, const vec3<uint32_t>* p_triangles, const mutable_adjacency_list* p_adj_list) :
    trivial_ear(p0, p1, p2, p_vertices, p_triangles, p_adj_list)
    {
    compute_quality(p_triangles);
    }

  JTKGINLINE void minimum_weight_ear::compute_quality(const vec3<uint32_t>* triangles)
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
    float dihedral_a = dihedral_angle(vertices[v0], vertices[v1], vertices[v2], vertices[a]);
    float dihedral_b = dihedral_angle(vertices[v1], vertices[v2], vertices[v0], vertices[b]);
    dihedral_rad = dihedral_a > dihedral_b ? dihedral_a : dihedral_b;
    }

  JTKGINLINE bool minimum_weight_ear::operator < (const minimum_weight_ear& other) const
    {
    if (is_concave() && !other.is_concave())
      return true;
    if (!is_concave() && other.is_concave())
      return false;
    return (aspect_ratio - (dihedral_rad / 3.141592653589793238462643383f) * 0.1f) < (other.aspect_ratio - (other.dihedral_rad / 3.141592653589793238462643383f) * 0.1f);
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

  } // namespace jtk

#endif // #ifndef JTK_GEOMETRY_H


#ifdef JTK_GEOMETRY_IMPLEMENTATION


#include "concurrency.h"
//#include "vec.h"
#include "point_tree.h"
#include <array>
#include <list>
#include <numeric>
#include <iterator>
#include <queue>
#include <set>
//#include <string>
//#include <unordered_map>
//#include <vector>

  /////////////////////////////////////////////////////////////////////////
  // implementation
  /////////////////////////////////////////////////////////////////////////


namespace jtk
  {

  JTKGDEF bool read_stl(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
    {
    FILE* inputfile;

    inputfile = fopen(filename, "rb");

    if (!inputfile)
      return false;

    char buffer[80];
    fread(buffer, 1, 80, inputfile);

    if (buffer[0] == 's' && buffer[1] == 'o' && buffer[2] == 'l' && buffer[3] == 'i' && buffer[4] == 'd')
      {
      // check whether it is an ascii file
      fclose(inputfile);
      std::ifstream inputfilecheck;
      inputfilecheck.open(filename);
      if (inputfilecheck.fail())
        return false;
      if (inputfilecheck.eof())
        return false;
      if (!inputfilecheck.is_open())
        return false;      
      std::string line1, line2;
      std::getline(inputfilecheck, line1);
      std::getline(inputfilecheck, line2);
      inputfilecheck.close();
      std::stringstream ss;
      ss << line2;
      std::string w1, w2;
      ss >> w1 >> w2;
      if (w1 == std::string("facet") && w2 == std::string("normal")) // this is an ascii file
        return false;
      // binary stl starting with solid in header, so continue reading
      inputfile = fopen(filename, "rb");
      if (!inputfile)
        return false;
      fread(buffer, 1, 80, inputfile);
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

    parallel_sort(vert.begin(), vert.end(), less_fie);

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

  JTKGDEF bool read_stl_ascii(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
    {
    std::ifstream inputfile;
    inputfile.open(filename);
    if (inputfile.fail())
      return false;
    if (inputfile.eof())
      return false;
    if (!inputfile.is_open())
      return false;

    {
    std::vector<vec3<float>>().swap(vertices);
    }
    {
    std::vector<vec3<uint32_t>>().swap(triangles);
    }

    std::vector<vec4<float>> vert;
    uint32_t nr_of_vertices = 0;

    while (!inputfile.eof())
      {
      std::string line;
      std::getline(inputfile, line);
      std::transform(line.begin(), line.end(), line.begin(), [](char ch) {return (char)::tolower(ch); });
      std::stringstream ss;
      ss << line;
      std::string w1, w2;
      ss >> w1 >> w2;
      if (w1 == std::string("outer") && w2 == std::string("loop"))
        {
        if (inputfile.eof())
          return false;
        triangles.emplace_back();
        std::getline(inputfile, line);
        float x, y, z;
        std::stringstream ss1;
        ss1 << line;
        ss1 >> w1 >> x >> y >> z;
        if (w1 != std::string("vertex"))
          return false;
        vert.emplace_back(x, y, z, *reinterpret_cast<const float*>(&nr_of_vertices));
        triangles.back()[0] = nr_of_vertices++;
        if (inputfile.eof())
          return false;
        std::getline(inputfile, line);
        std::stringstream ss2;
        ss2 << line;
        ss2 >> w1 >> x >> y >> z;
        if (w1 != std::string("vertex"))
          return false;
        vert.emplace_back(x, y, z, *reinterpret_cast<const float*>(&nr_of_vertices));
        triangles.back()[1] = nr_of_vertices++;
        if (inputfile.eof())
          return false;
        std::getline(inputfile, line);
        std::stringstream ss3;
        ss3 << line;
        ss3 >> w1 >> x >> y >> z;
        if (w1 != std::string("vertex"))
          return false;
        vert.emplace_back(x, y, z, *reinterpret_cast<const float*>(&nr_of_vertices));
        triangles.back()[2] = nr_of_vertices++;
        if (inputfile.eof())
          return false;
        std::getline(inputfile, line);
        std::stringstream ss4;
        ss4 << line;
        ss4 >> w1;
        if (w1 != std::string("endloop"))
          return false;
        }
      }

    inputfile.close();

    if (vert.empty())
      return false;

    auto less_fie = [](const vec4<float>& left, const vec4<float>& right)
      {
      return left[0] == right[0] ? (left[1] == right[1] ? left[2] < right[2] : left[1] < right[1]) : left[0] < right[0];
      };

    auto equal_fie = [](const vec4<float>& left, const vec4<float>& right)
      {
      return (left[0] == right[0] && left[1] == right[1] && left[2] == right[2]);
      };

    parallel_sort(vert.begin(), vert.end(), less_fie);

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

  JTKGDEF bool write_stl(const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const vec3<float>* triangle_normals, const unsigned short* attributes, const char* filename)
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

  JTKGDEF bool read_off(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
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

  JTKGDEF bool write_off(uint32_t nr_of_vertices, const vec3<float>* vertices, uint32_t nr_of_triangles, const vec3<uint32_t>* triangles, const char* filename)
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

  JTKGDEF bool read_obj(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, const char* filename)
    {
    FILE* f = nullptr;
    f = fopen(filename, "r");
    if (!f)
      return false;


    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      int first_non_whitespace_index = 0;
      while (first_non_whitespace_index < 250 && (buffer[first_non_whitespace_index] == ' ' || buffer[first_non_whitespace_index] == '\t'))
        ++first_non_whitespace_index;

      if (buffer[first_non_whitespace_index + 0] == 'v' && buffer[first_non_whitespace_index + 1] == ' ')
        {
        float x, y, z;
        auto err = sscanf(buffer + first_non_whitespace_index, "v %f %f %f\n", &x, &y, &z);
        if (err != 3)
          {
          fclose(f);
          return false;
          }
        vertices.push_back(vec3<float>(x, y, z));
        }
      else if (buffer[first_non_whitespace_index + 0] == 'f' && buffer[first_non_whitespace_index + 1] == ' ')
        {
        uint32_t t0, t1, t2, t3, v0, v1, v2, v3;
        auto err = sscanf(buffer + first_non_whitespace_index, "f %d/%d %d/%d %d/%d %d/%d\n", &t0, &v0, &t1, &v1, &t2, &v2, &t3, &v3);
        if (err == 6)
          {
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          }
        else if (err == 8)
          {
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
          }
        else
          {
          err = sscanf(buffer + first_non_whitespace_index, "f %d//%d %d//%d %d//%d %d//%d\n", &t0, &v0, &t1, &v1, &t2, &v2, &t3, &v3);
          if (err == 6)
            {
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            }
          else if (err == 8)
            {
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
            }
          else
            {
            err = sscanf(buffer + first_non_whitespace_index, "f %d %d %d %d\n", &t0, &t1, &t2, &t3);
            if (err == 3)
              {
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              }
            else if (err == 4)
              {
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
              }
            else
              {
              uint32_t tx0, tx1, tx2, tx3;
              err = sscanf(buffer + first_non_whitespace_index, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n", &t0, &v0, &tx0, &t1, &v1, &tx1, &t2, &v2, &tx2, &t3, &v3, &tx3);
              if (err == 9)
                {
                triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
                }
              else if (err == 12)
                {
                triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
                triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
                }
              else
                {
                fclose(f);
                return false;
                }

              }
            }
          }
        }
      }
    fclose(f);
    return true;
    }

  JTKGDEF bool read_texture_filename_from_mtl(std::string& texture_file, const char* filename)
    {
    FILE* f = nullptr;
    f = fopen(filename, "r");
    if (!f)
      return false;
    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      int first_non_whitespace_index = 0;
      while (first_non_whitespace_index < 250 && (buffer[first_non_whitespace_index] == ' ' || buffer[first_non_whitespace_index] == '\t'))
        ++first_non_whitespace_index;
      if (buffer[first_non_whitespace_index + 0] == 'm' && buffer[first_non_whitespace_index + 1] == 'a' && buffer[first_non_whitespace_index + 2] == 'p' && buffer[first_non_whitespace_index + 3] == '_' && buffer[first_non_whitespace_index + 4] == 'K' && buffer[first_non_whitespace_index + 5] == 'd')
        {
        char texture[256];
        auto scan_err = sscanf(buffer + first_non_whitespace_index, "map_Kd %s\n", texture);
        (void)scan_err;
        texture_file = std::string(texture);
        fclose(f);
        return true;
        }
      }
    fclose(f);
    return false;
    }

  JTKGDEF bool read_obj(std::string& mtl_filename, std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, std::vector<vec3<vec2<float>>>& uv, const char* filename)
    {
    mtl_filename = "";
    std::vector<vec3<float>>().swap(vertices);
    std::vector<vec3<uint32_t>>().swap(triangles);
    std::vector<vec3<vec2<float>>>().swap(uv);
    FILE* f = nullptr;
    f = fopen(filename, "r");
    if (!f)
      return false;

    std::vector<vec2<float>> tex;
    std::vector<vec3<uint32_t>> tria_uv;

    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      int first_non_whitespace_index = 0;
      while (first_non_whitespace_index < 250 && (buffer[first_non_whitespace_index] == ' ' || buffer[first_non_whitespace_index] == '\t'))
        ++first_non_whitespace_index;

      if (buffer[first_non_whitespace_index + 0] == 'm' && buffer[first_non_whitespace_index + 1] == 't' && buffer[first_non_whitespace_index + 2] == 'l' && buffer[first_non_whitespace_index + 3] == 'l' && buffer[first_non_whitespace_index + 4] == 'i' && buffer[first_non_whitespace_index + 5] == 'b')
        {
        char mtlfile[256];
        auto scan_err = sscanf(buffer + first_non_whitespace_index, "mtllib %s\n", mtlfile);
        if (scan_err == 1)
          {
          mtl_filename = std::string(mtlfile);
          }
        }
      if (buffer[first_non_whitespace_index + 0] == 'v' && buffer[first_non_whitespace_index + 1] == ' ')
        {
        float x, y, z;
        auto err = sscanf(buffer + first_non_whitespace_index, "v %f %f %f\n", &x, &y, &z);
        if (err != 3)
          {
          fclose(f);
          return false;
          }
        vertices.push_back(vec3<float>(x, y, z));
        }
      else if (buffer[first_non_whitespace_index + 0] == 'v' && buffer[first_non_whitespace_index + 1] == 't' && buffer[2] == ' ')
        {
        float x, y;
        auto err = sscanf(buffer + first_non_whitespace_index, "vt %f %f\n", &x, &y);
        if (err != 2)
          {
          fclose(f);
          return false;
          }
        tex.push_back(vec2<float>(x, y));
        }
      else if (buffer[first_non_whitespace_index + 0] == 'f' && buffer[first_non_whitespace_index + 1] == ' ')
        {
        uint32_t t0, t1, t2, t3, v0, v1, v2, v3;
        auto err = sscanf(buffer + first_non_whitespace_index, "f %d/%d %d/%d %d/%d %d/%d\n", &t0, &v0, &t1, &v1, &t2, &v2, &t3, &v3);
        if (err == 6)
          {
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          }
        else if (err == 8)
          {
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v2 - 1, v3 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
          }
        else
          {
          err = sscanf(buffer + first_non_whitespace_index, "f %d//%d %d//%d %d//%d %d//%d\n", &t0, &v0, &t1, &v1, &t2, &v2, &t3, &v3);
          if (err == 6)
            {
            tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            }
          else if (err == 8)
            {
            tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            tria_uv.push_back(vec3<uint32_t>(v0 - 1, v2 - 1, v3 - 1));
            triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
            }
          else
            {
            err = sscanf(buffer + first_non_whitespace_index, "f %d %d %d %d\n", &t0, &t1, &t2, &t3);
            if (err == 3)
              {
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              }
            else if (err == 4)
              {
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
              }
            else
              {
              uint32_t tx0, tx1, tx2, tx3;
              err = sscanf(buffer + first_non_whitespace_index, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n", &t0, &v0, &tx0, &t1, &v1, &tx1, &t2, &v2, &tx2, &t3, &v3, &tx3);
              if (err == 9)
                {
                tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
                triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
                }
              else if (err == 12)
                {
                tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
                triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
                tria_uv.push_back(vec3<uint32_t>(v0 - 1, v2 - 1, v3 - 1));
                triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
                }
              else
                {
                fclose(f);
                return false;
                }

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

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts)
    {
    FILE* fp = fopen(filename, "wt");
    if (!fp)
      return false;
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
    return true;
    }

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<float>>& normals)
    {
    if (normals.empty())
      {
      return write_ply(filename, pts);
      }
    FILE* fp = fopen(filename, "wt");
    if (!fp)
      return false;
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
    return true;
    }

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<uint32_t>& clrs)
    {
    if (clrs.empty())
      {
      return write_ply(filename, pts);
      }
    FILE* fp = fopen(filename, "wt");
    if (!fp)
      return false;
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
    return true;
    }

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<float>>& normals, const std::vector<uint32_t>& clrs)
    {
    if (normals.empty())
      {
      return write_ply(filename, pts, clrs);
      }
    if (clrs.empty())
      {
      return write_ply(filename, pts, normals);
      }
    FILE* fp = fopen(filename, "wt");
    if (!fp)
      return false;
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
    return true;
    }

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<uint32_t>>& triangles)
    {
    FILE* fp = fopen(filename, "wb");
    if (!fp)
      return false;
    fprintf(fp, "ply\n");
    int n = 1;
    if (*(char*)&n == 1)
      fprintf(fp, "format binary_little_endian 1.0\n");
    else
      fprintf(fp, "format binary_big_endian 1.0\n");
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
    return true;
    }


  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<uint32_t>& pts_colors, const std::vector<vec3<uint32_t>>& triangles)
    {
    FILE* fp = fopen(filename, "wb");
    if (!fp)
      return false;
    fprintf(fp, "ply\n");
    int n = 1;
    if (*(char*)&n == 1)
      fprintf(fp, "format binary_little_endian 1.0\n");
    else
      fprintf(fp, "format binary_big_endian 1.0\n");
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
    return true;
    }

  JTKGDEF bool write_ply(const char* filename, const std::vector<vec3<float>>& pts, const std::vector<vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<jtk::vec2<float>>>& uv)
    {
    FILE* fp = fopen(filename, "wb");
    if (!fp)
      return false;
    fprintf(fp, "ply\n");
    int n = 1;
    if (*(char*)&n == 1)
      fprintf(fp, "format binary_little_endian 1.0\n");
    else
      fprintf(fp, "format binary_big_endian 1.0\n");
    fprintf(fp, "element vertex %d\n", (int)pts.size());
    fprintf(fp, "property float x\n");
    fprintf(fp, "property float y\n");
    fprintf(fp, "property float z\n");
    fprintf(fp, "element face %d\n", (int)triangles.size());
    fprintf(fp, "property list uchar int vertex_indices\n");
    fprintf(fp, "property list uchar float texcoord\n");
    fprintf(fp, "end_header\n");
    fwrite(pts.data(), sizeof(vec3<float>), pts.size(), fp);
    unsigned char tria_size = 3;
    const unsigned char texcoord_size = 6;
    for (uint32_t i = 0; i < (uint32_t)triangles.size(); ++i)
      {
      fwrite(&tria_size, 1, 1, fp);
      fwrite((uint32_t*)triangles.data() + 3 * i, sizeof(uint32_t), 3, fp);
      if (!uv.empty())
        {
        fwrite(&texcoord_size, 1, 1, fp);
        fwrite((float*)uv.data() + 6 * i, sizeof(float), 6, fp);
        }
      }
    fclose(fp);
    return true;
    }

  JTKGDEF void delete_triangles(std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& triangles_to_delete)
    {
    delete_items(triangles, triangles_to_delete);
    }

  JTKGDEF void remove_free_vertices(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles)
    {
    std::vector<uint32_t> count_occurrence(vertices.size(), 0);
    for (const auto& tria : triangles)
      {
      ++count_occurrence[tria[0]];
      ++count_occurrence[tria[1]];
      ++count_occurrence[tria[2]];
      }
    std::vector<uint32_t> diff(vertices.size(), 0);
    uint32_t current_off = 0;
    std::vector<vec3<float>> new_vertices;
    auto cnt = (size_t)std::count(count_occurrence.begin(), count_occurrence.end(), (uint32_t)0);
    if (cnt <= vertices.size())
      {
      new_vertices.reserve((size_t)(vertices.size() - cnt));
      }
    for (size_t v = 0; v < vertices.size(); ++v)
      {
      if (count_occurrence[v] == 0)
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

  JTKGDEF void remove_free_vertices(std::vector<vec3<float>>& vertices, std::vector<uint32_t>& colors, std::vector<vec3<uint32_t>>& triangles)
    {
    std::vector<uint32_t> count_occurrence(vertices.size(), 0);
    for (const auto& tria : triangles)
      {
      ++count_occurrence[tria[0]];
      ++count_occurrence[tria[1]];
      ++count_occurrence[tria[2]];
      }
    std::vector<uint32_t> diff(vertices.size(), 0);
    uint32_t current_off = 0;
    std::vector<vec3<float>> new_vertices;
    std::vector<uint32_t> new_colors;
    auto cnt = (size_t)std::count(count_occurrence.begin(), count_occurrence.end(), (uint32_t)0);
    if (cnt <= vertices.size())
      {
      new_vertices.reserve((size_t)(vertices.size() - cnt));
      new_colors.reserve((size_t)(vertices.size() - cnt));
      }
    for (size_t v = 0; v < vertices.size(); ++v)
      {
      if (count_occurrence[v] == 0)
        {
        ++current_off;
        }
      else
        {
        new_vertices.push_back(vertices[v]);
        new_colors.push_back(colors[v]);
        }
      diff[v] = current_off;
      }

    for (auto& tria : triangles)
      {
      tria[0] -= diff[tria[0]];
      tria[1] -= diff[tria[1]];
      tria[2] -= diff[tria[2]];
      }
    vertices.swap(new_vertices);
    colors.swap(new_colors);
    }

  JTKGDEF bool edge_swap(uint32_t v, uint32_t v2, std::vector<jtk::vec3<uint32_t>>& triangles, jtk::mutable_adjacency_list& adj_list)
    {
    auto tria = jtk::triangle_indices_from_edge(v, v2, adj_list);
    if (tria.size() != 2)
      return false;

    auto t1 = triangles[tria[0]];
    auto t2 = triangles[tria[1]];
    int t1v = 0;
    while (t1[t1v] == v || t1[t1v] == v2)
      ++t1v;
    int t2v = 0;
    while (t2[t2v] == v || t2[t2v] == v2)
      ++t2v;

    int vpos = 0;
    while (t1[vpos] != v)
      ++vpos;
    int v2pos = 0;
    while (t1[v2pos] != v2)
      ++v2pos;

    if ((v2pos + 1) % 3 == vpos)
      {
      std::swap(tria[0], tria[1]);
      std::swap(t1, t2);
      std::swap(t1v, t2v);
      }

    uint32_t T1v = t1[t1v];
    uint32_t T2v = t2[t2v];

    if (!jtk::triangle_indices_from_edge(T1v, T2v, adj_list).empty())
      return false;

    adj_list.remove_triangle_from_vertex(v2, tria[0]);
    adj_list.remove_triangle_from_vertex(v, tria[1]);

    t1[0] = T1v;
    t1[1] = v;
    t1[2] = T2v;

    t2[0] = T1v;
    t2[1] = T2v;
    t2[2] = v2;

    triangles[tria[0]] = t1;
    triangles[tria[1]] = t2;

    adj_list.add_triangle_to_vertex(T1v, tria[1]);
    adj_list.add_triangle_to_vertex(T2v, tria[0]);

    assert(jtk::triangle_indices_from_edge(T1v, T2v, adj_list).size() == 2);
    assert(jtk::triangle_indices_from_edge(T1v, v2, adj_list).size() <= 2);
    assert(jtk::triangle_indices_from_edge(v2, T2v, adj_list).size() <= 2);
    assert(jtk::triangle_indices_from_edge(T1v, v, adj_list).size() <= 2);
    assert(jtk::triangle_indices_from_edge(v, T2v, adj_list).size() <= 2);
    return true;
    }

  JTKGDEF std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles)
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

  JTKGDEF std::vector<uint32_t> one_ring_vertices_from_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles)
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

  JTKGDEF std::vector<std::vector<uint32_t>> ordered_one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles, bool oriented)
    {
    auto it = adj_list.begin(vertex_index);
    const auto it_end = adj_list.end(vertex_index);
    if (it == it_end)
      return std::vector<std::vector<uint32_t>>();

    std::list<uint32_t> trias;
    for (auto it2 = it; it2 != it_end; ++it2)
      {
      trias.push_back(*it2);
      }

    std::vector<std::vector<uint32_t>> batches;

    while (!trias.empty())
      {
      std::vector<uint32_t> current_loop;
      const auto tria = triangles[trias.front()];
      trias.pop_front();
      uint32_t v = 0;
      if (tria[v] != vertex_index)
        {
        ++v;
        if (tria[v] != vertex_index)
          ++v;
        }
      current_loop.push_back(tria[(v + 1) % 3]);
      current_loop.push_back(tria[(v + 2) % 3]);
      bool done = false;
      while (!done)
        {
        done = true;
        for (auto it3 = trias.begin(); it3 != trias.end(); ++it3)
          {
          const auto local_tria = triangles[*it3];
          v = 0;
          if (local_tria[v] != vertex_index)
            {
            ++v;
            if (local_tria[v] != vertex_index)
              ++v;
            }
          if (local_tria[(v + 1) % 3] == current_loop.back())
            {
            done = false;
            current_loop.push_back(local_tria[(v + 2) % 3]);
            trias.erase(it3);
            break;
            }
          else if (local_tria[(v + 2) % 3] == current_loop.front())
            {
            done = false;
            current_loop.insert(current_loop.begin(), local_tria[(v + 1) % 3]);
            trias.erase(it3);
            break;
            }
          else if (!oriented && (local_tria[(v + 2) % 3] == current_loop.back()))
            {
            done = false;
            current_loop.push_back(local_tria[(v + 1) % 3]);
            trias.erase(it3);
            break;
            }
          else if (!oriented && (local_tria[(v + 1) % 3] == current_loop.front()))
            {
            done = false;
            current_loop.insert(current_loop.begin(), local_tria[(v + 2) % 3]);
            trias.erase(it3);
            break;
            }

          }
        }
      if (current_loop.back() == current_loop.front())
        current_loop.pop_back();
      batches.push_back(current_loop);
      }
    return batches;
    }

  JTKGDEF vec3<float> normal(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles)
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

  JTKGDEF vec3<float> normal(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle)
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

  JTKGDEF void compute_triangle_normals(std::vector<vec3<float>>& triangle_normals, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles)
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

  JTKGDEF void compute_vertex_normals(std::vector<jtk::vec3<float>>& vertex_normals, const jtk::vec3<float>* triangle_normals,
    const jtk::vec3<float>* vertices, const uint32_t nr_of_vertices, const jtk::vec3<uint32_t>* triangles, const uint32_t nr_of_triangles)
    {
    vertex_normals.resize(nr_of_vertices, vec3<float>(0.f, 0.f, 0.f));

    std::vector<uint32_t> vertex_occurrence(nr_of_vertices + 1, 0);
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      ++vertex_occurrence[triangles[t][0]];
      ++vertex_occurrence[triangles[t][1]];
      ++vertex_occurrence[triangles[t][2]];
      }

    uint32_t prev_value = vertex_occurrence[0];
    vertex_occurrence[0] = 0;
    for (uint32_t v = 0; v < nr_of_vertices; ++v)
      {
      uint32_t current = vertex_occurrence[v + 1];
      vertex_occurrence[v + 1] = vertex_occurrence[v] + prev_value;
      prev_value = current;
      }

    std::vector<uint32_t> vertex_count(nr_of_vertices, 0);
    std::vector<uint32_t> triangle_list(vertex_occurrence.back());
    for (uint32_t t = 0; t < nr_of_triangles; ++t)
      {
      for (int j = 0; j < 3; ++j)
        {
        const uint32_t v = triangles[t][j];
        const uint32_t pos = vertex_occurrence[v];
        const uint32_t off = vertex_count[v]++;
        triangle_list[pos + off] = t;
        }
      }

    for (uint32_t v = 0; v < nr_of_vertices; ++v)
      {
      uint32_t* it = triangle_list.data() + vertex_occurrence[v];
      const uint32_t* it_end = triangle_list.data() + vertex_occurrence[v + 1];
      auto vn = vec3<float>(0.f);
      for (; it != it_end; ++it)
        {
        const uint32_t t = *it;
        const auto nn = triangle_normals[t];
        auto v0 = triangles[t][0];
        auto v1 = triangles[t][1];
        auto v2 = triangles[t][2];
        vec3<float> V0(vertices[v0][0], vertices[v0][1], vertices[v0][2]);
        vec3<float> V1(vertices[v1][0], vertices[v1][1], vertices[v1][2]);
        vec3<float> V2(vertices[v2][0], vertices[v2][1], vertices[v2][2]);
        float twice_triangle_length = length(jtk::cross(V1 - V0, V2 - V0));
        vn = vn + nn * twice_triangle_length;
        }
      auto length = std::sqrt(dot(vn, vn));
      if (length)
        vn = vn / length;
      vertex_normals[v] = vn;
      }
    }

  JTKGDEF std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list)
    {
    assert(v0 != v1);
    std::vector<uint32_t> result;
    std::set_intersection(adj_list.begin(v0), adj_list.end(v0), adj_list.begin(v1), adj_list.end(v1), std::back_inserter(result));
    return result;
    }

  JTKGDEF std::vector<uint32_t> triangle_indices_from_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list)
    {
    assert(v0 != v1);
    std::vector<uint32_t> result;
    result.reserve(2);
    std::set_intersection(adj_list.begin(v0), adj_list.end(v0), adj_list.begin(v1), adj_list.end(v1), std::back_inserter(result));
    return result;
    }

  JTKGDEF float signed_volume(const vec3<float>* vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles)
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

  JTKGDEF float area(uint32_t triangle_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles)
    {
    const auto v0 = triangles[triangle_index][0];
    const auto v1 = triangles[triangle_index][1];
    const auto v2 = triangles[triangle_index][2];
    const vec3<float> V0(vertices[v0]);
    const vec3<float> V1(vertices[v1]);
    const vec3<float> V2(vertices[v2]);
    return length(cross(V1 - V0, V2 - V0)) / 2.f;
    }

  JTKGDEF bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle)
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

  JTKGDEF bool is_sharp_edge(uint32_t v0, uint32_t v1, const vec3<float>* triangle_normals, const adjacency_list& adj_list, float cos_sharp_angle)
    {
    auto trias = triangle_indices_from_edge(v0, v1, adj_list);
    if (trias.size() == 2)
      {
      auto ca = dot(triangle_normals[trias[0]], triangle_normals[trias[1]]);
      return ca < cos_sharp_angle;
      }
    return false;
    }

  JTKGDEF bool vertex_on_sharp_edge(uint32_t vertex_index, const vec3<float>* vertices, const vec3<uint32_t>* triangles, const adjacency_list& adj_list, float cos_sharp_angle)
    {
    auto one_ring = one_ring_vertices_from_vertex(vertex_index, adj_list, triangles);
    for (auto v : one_ring)
      {
      if (is_sharp_edge(vertex_index, v, vertices, triangles, adj_list, cos_sharp_angle))
        return true;
      }
    return false;
    }

  JTKGDEF std::vector<std::vector<uint32_t>> triangle_neighbour_indices(const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, const adjacency_list& adj_list)
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

  JTKGDEF bool edge_between_triangles(uint32_t& v0, uint32_t& v1, const vec3<uint32_t>& tria1, const vec3<uint32_t>& tria2)
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

  JTKGDEF bool is_manifold_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list)
    {
    auto tri_idces = triangle_indices_from_edge(v0, v1, adj_list);
    return tri_idces.size() == 2;
    }

  JTKGDEF bool is_boundary_edge(uint32_t v0, uint32_t v1, const adjacency_list& adj_list)
    {
    return triangle_indices_from_edge(v0, v1, adj_list).size() == 1;
    }

  JTKGDEF bool is_boundary_edge(uint32_t v0, uint32_t v1, const mutable_adjacency_list& adj_list)
    {
    return triangle_indices_from_edge(v0, v1, adj_list).size() == 1;
    }

  JTKGDEF bool is_boundary_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    auto vert = one_ring_vertices_from_vertex(vertex_index, adj_list, triangles);
    for (auto v : vert)
      {
      if (is_boundary_edge(v, vertex_index, adj_list))
        return true;
      }
    return false;
    }

  JTKGDEF bool is_boundary_vertex(uint32_t vertex_index, const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles)
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
    JTKGDEF std::vector<uint32_t> _shells(const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, bool manifold)
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


  JTKGDEF std::vector<uint32_t> shells(const uint32_t nr_of_vertices, const vec3<uint32_t>* triangles, const uint32_t nr_of_triangles, shell_connectivity conn)
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

  JTKGDEF float dihedral_angle(const vec3<float>& u, const vec3<float>& v, const vec3<float>& a, const vec3<float>& b)
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

  namespace hole_filling_details
    {
    template <class TAdjList>
    std::vector<std::vector<uint32_t>> _holes(const TAdjList& adj_list, const vec3<uint32_t>* triangles)
      {
      std::vector<bool> vertex_treated(adj_list.size(), false);
      std::vector<std::vector<uint32_t>> holes;
      for (uint32_t v = 0; v < (uint32_t)vertex_treated.size(); ++v)
        {
        if (!vertex_treated[v])
          {
          vertex_treated[v] = true;
          if (is_boundary_vertex(v, adj_list, triangles))
            {
            std::vector<uint32_t> hole;
            hole.push_back(v);
            std::queue<uint32_t> qu;
            auto neighbouring_vertices = one_ring_vertices_from_vertex(v, adj_list, triangles);
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
              neighbouring_vertices = one_ring_vertices_from_vertex(current_vertex, adj_list, triangles);
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
              neighbouring_vertices = one_ring_vertices_from_vertex(v, adj_list, triangles);
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
      return holes;
      }
    }

  JTKGDEF std::vector<std::vector<uint32_t>> holes(const adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    return hole_filling_details::_holes<adjacency_list>(adj_list, triangles);
    }

  JTKGDEF std::vector<std::vector<uint32_t>> holes(const mutable_adjacency_list& adj_list, const vec3<uint32_t>* triangles)
    {
    return hole_filling_details::_holes<mutable_adjacency_list>(adj_list, triangles);
    }

  template <class Ear>
  void fill_holes(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, uint32_t max_holes_size)
    {
    mutable_adjacency_list adj_list((uint32_t)vertices.size(), triangles.data(), (uint32_t)triangles.size());

    std::vector<std::vector<uint32_t>> holes_found = holes(adj_list, triangles.data());

    for (auto& hole : holes_found)
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

  template void fill_hole_ear<minimum_weight_ear>(std::vector<vec3<uint32_t>>&, const vec3<float>*, mutable_adjacency_list&, const std::vector<uint32_t>&);
  template void fill_holes<minimum_weight_ear>(std::vector<vec3<float>>&, std::vector<vec3<uint32_t>>&, uint32_t);
  template void fill_hole_ear<trivial_ear>(std::vector<vec3<uint32_t>>&, const vec3<float>*, mutable_adjacency_list&, const std::vector<uint32_t>&);
  template void fill_holes<trivial_ear>(std::vector<vec3<float>>&, std::vector<vec3<uint32_t>>&, uint32_t);

  namespace details
    {

    JTKGDEF void _CreateFilter(std::vector<std::vector<uint32_t>>& io_filter, const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles)
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

    JTKGDEF void _ApplyFilter(std::vector<vec3<float>>& o_vertices, std::vector<vec3<float>>& i_vertices, const std::vector<std::vector<uint32_t>>& i_filter, float filter_value)
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
            o_vertices[j][k] += filter_value / float(valence) * i_vertices[i_filter[j][v]][k];
          }
        });
      }

    JTKGDEF void _ApplyLocalFilter(std::vector<vec3<float>>& o_vertices, std::vector<vec3<float>>& i_vertices, const std::vector<uint32_t>& vertex_indices, const std::vector<std::vector<uint32_t>>& i_filter, float filter_value)
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
              o_vertices[vertex_index][k] += filter_value / float(valence) * i_vertices[i_filter[vertex_index][v]][k];
            }
        });
      }
    }

  JTKGDEF void smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, uint32_t iterations, float lambda, float mu)
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

  JTKGDEF void local_smooth(std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, const std::vector<uint32_t>& vertex_indices, uint32_t iterations, float lambda, float mu)
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

  JTKGDEF void dyadic_subdivide(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles)
    {
    const uint32_t nr_of_vertices = (uint32_t)vertices.size();
    const uint32_t nr_of_triangles = (uint32_t)triangles.size();
    vertices.reserve(nr_of_vertices * 2 + nr_of_triangles + 10);

    std::unordered_map<uint64_t, uint32_t> treated(nr_of_vertices);

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

  JTKGDEF void dyadic_subdivide_uv_map(std::vector<jtk::vec3<jtk::vec2<float>>>& triangle_uv)
    {
    size_t original_triangles = triangle_uv.size();
    std::vector<jtk::vec3<jtk::vec2<float>>> new_uv_out;
    new_uv_out.resize(original_triangles * 4);
    for (size_t t = 0; t < original_triangles; ++t)
      {
      auto uv = triangle_uv[t];
      jtk::vec3<jtk::vec2<float>> new_uv;
      new_uv[0] = (uv[0] + uv[1]) * 0.5f;
      new_uv[1] = (uv[1] + uv[2]) * 0.5f;
      new_uv[2] = (uv[2] + uv[0]) * 0.5f;
      new_uv_out[t] = new_uv;
      new_uv[0] = uv[0];
      new_uv[1] = (uv[0] + uv[1]) * 0.5f;
      new_uv[2] = (uv[2] + uv[0]) * 0.5f;
      new_uv_out[original_triangles + 3 * t] = new_uv;
      new_uv[1] = uv[1];
      new_uv[2] = (uv[1] + uv[2]) * 0.5f;
      new_uv[0] = (uv[0] + uv[1]) * 0.5f;
      new_uv_out[original_triangles + 3 * t + 1] = new_uv;
      new_uv[2] = uv[2];
      new_uv[0] = (uv[2] + uv[0]) * 0.5f;
      new_uv[1] = (uv[1] + uv[2]) * 0.5f;
      new_uv_out[original_triangles + 3 * t + 2] = new_uv;
      }
    triangle_uv.swap(new_uv_out);
    }

  JTKGDEF void dyadic_subdivide_triangle_indices_vector(std::vector<uint32_t>& triangle_indices, uint32_t original_number_of_triangles)
    {
    const size_t sz = triangle_indices.size();
    triangle_indices.reserve(sz * 4);
    for (size_t i = 0; i < sz; ++i)
      {
      uint32_t t = triangle_indices[i];
      triangle_indices.push_back(original_number_of_triangles + t * 3);
      triangle_indices.push_back(original_number_of_triangles + t * 3 + 1);
      triangle_indices.push_back(original_number_of_triangles + t * 3 + 2);
      }
    }

  JTKGDEF void undo_dyadic_subdivide(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles)
    {
    const uint32_t nr_of_vertices = (uint32_t)vertices.size();
    const uint32_t nr_of_triangles = (uint32_t)triangles.size();

    assert(nr_of_triangles % 4 == 0);

    uint32_t original_nr_of_vertices = triangles.front()[0];
    const uint32_t original_nr_of_triangles = nr_of_triangles / 4;
    for (uint32_t t = 0; t < original_nr_of_triangles; ++t)
      {
      original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][0]);
      original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][1]);
      original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][2]);
      }
    for (uint32_t t = 0; t < original_nr_of_triangles; ++t)
      {
      const uint32_t t0 = original_nr_of_triangles + t * 3;
      const uint32_t t1 = original_nr_of_triangles + t * 3 + 1;
      const uint32_t t2 = original_nr_of_triangles + t * 3 + 2;
      triangles[t][0] = triangles[t0][0];
      triangles[t][1] = triangles[t1][1];
      triangles[t][2] = triangles[t2][2];
      }
    triangles.resize(original_nr_of_triangles);
    vertices.resize(original_nr_of_vertices);
    }

  JTKGDEF void undo_dyadic_subdivide_uv_map(std::vector<jtk::vec3<jtk::vec2<float>>>& triangle_uv)
    {
    const uint32_t number_of_original_triangles = (uint32_t)triangle_uv.size() / 4;
    for (uint32_t t = 0; t < number_of_original_triangles; ++t)
      {
      triangle_uv[t][0] = triangle_uv[number_of_original_triangles + 3 * t][0];
      triangle_uv[t][1] = triangle_uv[number_of_original_triangles + 3 * t + 1][1];
      triangle_uv[t][2] = triangle_uv[number_of_original_triangles + 3 * t + 2][2];
      }
    triangle_uv.resize(number_of_original_triangles);
    }

  JTKGDEF void butterfly(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles)
    {
    const uint32_t nr_of_vertices = (uint32_t)vertices.size();
    const uint32_t nr_of_triangles = (uint32_t)triangles.size();

    assert(nr_of_triangles % 4 == 0);

    uint32_t original_nr_of_vertices = triangles.front()[0];
    const uint32_t original_nr_of_triangles = nr_of_triangles / 4;
    for (uint32_t t = 0; t < original_nr_of_triangles; ++t)
      {
      original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][0]);
      original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][1]);
      original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][2]);
      }

    std::vector<jtk::vec3<uint32_t>> old_triangles;
    old_triangles.reserve(original_nr_of_triangles);
    for (uint32_t t = 0; t < original_nr_of_triangles; ++t)
      {
      const uint32_t t0 = original_nr_of_triangles + t * 3;
      const uint32_t t1 = original_nr_of_triangles + t * 3 + 1;
      const uint32_t t2 = original_nr_of_triangles + t * 3 + 2;
      old_triangles.emplace_back(triangles[t0][0], triangles[t1][1], triangles[t2][2]);
      }

    adjacency_list new_adj_list(nr_of_vertices, triangles.data(), nr_of_triangles);
    adjacency_list old_adj_list(original_nr_of_vertices, old_triangles.data(), original_nr_of_triangles);
    for (uint32_t v = original_nr_of_vertices; v < (uint32_t)vertices.size(); ++v)
      {
      auto neighbours = one_ring_vertices_from_vertex(v, new_adj_list, triangles.data());
      uint32_t e[2];
      int e_index = 0;
      for (const auto& n : neighbours)
        {
        if (n < original_nr_of_vertices)
          e[e_index++] = n;
        }

      if (is_boundary_edge(e[0], e[1], old_adj_list))
        {
        uint32_t e2[2];
        for (int k = 0; k < 2; ++k)
          {
          auto neighbours2 = one_ring_vertices_from_vertex(e[k], old_adj_list, old_triangles.data());
          for (const auto& n : neighbours2)
            {
            if (n != e[(k + 1) % 2] && is_boundary_edge(e[k], n, old_adj_list))
              {
              e2[k] = n;
              break;
              }
            }
          }
        auto new_vertex_position = (vertices[e2[0]] + vertices[e2[1]]) / -16.f + 9.f * (vertices[e[0]] + vertices[e[1]]) / 16.f;
        vertices[v] = new_vertex_position;
        }
      else
        {
        auto n0 = ordered_one_ring_vertices_from_vertex(e[0], old_adj_list, old_triangles.data()).front();
        auto n1 = ordered_one_ring_vertices_from_vertex(e[1], old_adj_list, old_triangles.data()).front();

        size_t e1_position = std::distance(n0.begin(), std::find(n0.begin(), n0.end(), e[1]));
        size_t e0_position = std::distance(n1.begin(), std::find(n1.begin(), n1.end(), e[0]));

        assert(e1_position < n0.size());
        assert(e0_position < n1.size());

        if (n0.size() == 6 && n1.size() == 6) // regular butterfly
          {
          uint32_t b0_up = n0[(e1_position + 1) % n0.size()];
          uint32_t b0_down = n0[(e1_position + n0.size() - 1) % n0.size()];
          uint32_t c0_up = n0[(e1_position + 2) % n0.size()];
          uint32_t c0_down = n0[(e1_position + n0.size() - 2) % n0.size()];

          //uint32_t b1_down = n1[(e0_position+1)%n1.size()];
          //uint32_t b1_up = n1[(e0_position+n1.size()-1)%n1.size()];
          uint32_t c1_down = n1[(e0_position + 2) % n1.size()];
          uint32_t c1_up = n1[(e0_position + n0.size() - 2) % n1.size()];

          assert(b0_up == n1[(e0_position + n1.size() - 1) % n1.size()]);
          assert(b0_down == n1[(e0_position + 1) % n1.size()]);

          auto new_vertex_position = (vertices[e[0]] + vertices[e[1]]) * 0.5f + (vertices[b0_up] + vertices[b0_down]) / 8.f - (vertices[c0_up] + vertices[c0_down] + vertices[c1_up] + vertices[c1_down]) / 16.f;
          vertices[v] = new_vertex_position;
          }
        else if (n0.size() == 6 || n1.size() == 6)
          {
          if (n0.size() == 6)
            {
            std::swap(n0, n1);
            std::swap(e0_position, e1_position);
            }
          std::vector<float> weights;
          weights.reserve(n0.size());
          if (n0.size() == 3)
            {
            weights.push_back(5.f / 12.f);
            weights.push_back(-1.f / 12.f);
            weights.push_back(-1.f / 12.f);
            }
          else if (n0.size() == 4)
            {
            weights.push_back(3.f / 8.f);
            weights.push_back(0.f);
            weights.push_back(-1.f / 8.f);
            weights.push_back(0.f);
            }
          else
            {
            for (int j = 0; j < (int)n0.size(); ++j)
              {
              weights.push_back((float)((0.25 + std::cos(6.28318530718f * (float)j / (float)n0.size()) + 0.5 * std::cos(12.5663706144f * (float)j / (float)n0.size())) / (float)n0.size()));
              }
            }

          float q = 1.f - std::accumulate(weights.begin(), weights.end(), 0.f);
          jtk::vec3<float> new_vertex_position = q * vertices[n1[e0_position]];
          for (int j = 0; j < (int)n0.size(); ++j)
            {
            new_vertex_position = new_vertex_position + vertices[n0[(j + e1_position) % n0.size()]] * weights[j];
            }
          vertices[v] = new_vertex_position;
          }
        else
          {
          std::vector<float> weights0;
          weights0.reserve(n0.size());
          if (n0.size() == 3)
            {
            weights0.push_back(5.f / 12.f);
            weights0.push_back(-1.f / 12.f);
            weights0.push_back(-1.f / 12.f);
            }
          else if (n0.size() == 4)
            {
            weights0.push_back(3.f / 8.f);
            weights0.push_back(0.f);
            weights0.push_back(-1.f / 8.f);
            weights0.push_back(0.f);
            }
          else
            {
            for (int j = 0; j < (int)n0.size(); ++j)
              {
              weights0.push_back((float)((0.25 + std::cos(6.28318530718f * (float)j / (float)n0.size()) + 0.5 * std::cos(12.5663706144f * (float)j / (float)n0.size())) / (float)n0.size()));
              }
            }

          float q0 = 1.f - std::accumulate(weights0.begin(), weights0.end(), 0.f);
          jtk::vec3<float> new_vertex_position_0 = q0 * vertices[n1[e0_position]];
          for (int j = 0; j < (int)n0.size(); ++j)
            {
            new_vertex_position_0 = new_vertex_position_0 + vertices[n0[(j + e1_position) % n0.size()]] * weights0[j];
            }

          std::vector<float> weights1;
          weights1.reserve(n1.size());
          if (n1.size() == 3)
            {
            weights1.push_back(5.f / 12.f);
            weights1.push_back(-1.f / 12.f);
            weights1.push_back(-1.f / 12.f);
            }
          else if (n1.size() == 4)
            {
            weights1.push_back(3.f / 8.f);
            weights1.push_back(0.f);
            weights1.push_back(-1.f / 8.f);
            weights1.push_back(0.f);
            }
          else
            {
            for (int j = 0; j < (int)n1.size(); ++j)
              {
              weights1.push_back((float)((0.25 + std::cos(6.28318530718f * (float)j / (float)n1.size()) + 0.5 * std::cos(12.5663706144f * (float)j / (float)n1.size())) / (float)n1.size()));
              }
            }

          float q1 = 1.f - std::accumulate(weights1.begin(), weights1.end(), 0.f);
          jtk::vec3<float> new_vertex_position_1 = q1 * vertices[n0[e1_position]];
          for (int j = 0; j < (int)n1.size(); ++j)
            {
            new_vertex_position_1 = new_vertex_position_1 + vertices[n1[(j + e0_position) % n1.size()]] * weights1[j];
            }

          vertices[v] = (new_vertex_position_0 + new_vertex_position_1) * 0.5f;
          }
        }
      }
    }

  JTKGDEF bool stitch_points(std::vector<vec3<float>>& vertices, std::vector<vec3<uint32_t>>& triangles, float distance_tolerance)
    {
    if (vertices.empty())
      return false;

    struct indexed_point
      {
      vec3<float> pt;
      uint32_t idx;
      float operator[](int i) const { return pt[i]; }
      float& operator[](int i) { return pt[i]; }
      };

    struct point_tree_traits
      {
      typedef float value_type;
      enum { dimension = 3 };
      typedef indexed_point point_type;
      };

    std::vector<indexed_point> points;
    points.resize(vertices.size());
    for (uint32_t i = 0; i < (uint32_t)vertices.size(); ++i)
      {
      points[i].pt = vertices[i];
      points[i].idx = i;
      }

    point_tree<point_tree_traits> tree;
    tree.efficient_build_tree(points.begin(), points.end());

    bool did_change = false;

    uint32_t number_of_vertices = (uint32_t)vertices.size();
    std::vector<uint32_t> new_indices(number_of_vertices, number_of_vertices);
    uint32_t new_index = 0;
    for (uint32_t i = 0; i < number_of_vertices; ++i)
      {
      if (new_indices[i] == number_of_vertices)
        {
        new_indices[i] = new_index;
        indexed_point pt;
        pt.pt = vertices[i];
        auto closest = tree.find_nearest_within_radius(distance_tolerance, pt);
        for (auto it = closest.begin(); it != closest.end(); ++it)
          if (new_indices[it->idx] == number_of_vertices && it->idx != i) {
            did_change = true;
            new_indices[it->idx] = new_index;
            }
        ++new_index;
        }
      }
    new_index = 0;
    for (uint32_t i = 0; i < number_of_vertices; ++i)
      {
      if (new_indices[i] == new_index)
        {
        vertices[new_index] = vertices[i];
        ++new_index;
        }
      }
    vertices.resize(new_index);

    for (auto& tria : triangles)
      for (uint32_t j = 0; j < 3; ++j)
        tria[j] = new_indices[tria[j]];

    return did_change;
    }

  } // namespace jtk

#endif// JTK_GEOMETRY_IMPLEMENTATION
