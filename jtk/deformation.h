/*
   Do this:
      #define JTK_DEFORMATION_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_DEFORMATION_IMPLEMENTATION
   #include "jtk/deformation.h"
 */

#ifndef JTK_DEFORMATION_H
#define JTK_DEFORMATION_H

#ifndef JTKDEFORMATIONDEF
#ifdef JTK_DEFORMATION_STATIC
#define JTK_QBVH_STATIC
#define JTK_GEOMETRY_STATIC
#define JTKDEFORMATIONDEF static
#else
#define JTKDEFORMATIONDEF extern
#endif
#endif

#ifndef JTKDEFORMATIONINLINE
#define JTKDEFORMATIONINLINE inline 
#endif

#include "vec.h"
#include "qbvh.h"
#include "geometry.h"
#include "concurrency.h"
#include <memory>



namespace jtk
  {

  template <class TDeformationTool>
  void deform(std::vector<vec3<float>>& pts, TDeformationTool& tool)
    {
    size_t sz = pts.size();
    parallel_for(size_t(0), sz, [&](size_t i)
      {
      if (tool.influences_point(pts[i]))
        tool.transform(pts[i]);
      });
    }

  template <class TDeformationTool>
  void deform(std::vector<vec3<float>>& pts, const float4x4& pts_coordinate_system, TDeformationTool& tool)
    {
    auto cs_inv = invert_orthonormal(pts_coordinate_system);
    size_t sz = pts.size();
    parallel_for(size_t(0), sz, [&](size_t i)
      {
      auto pt = transform(pts_coordinate_system, pts[i]);
      if (tool.influences_point(pt))
        tool.transform(pt);
      pts[i] = transform(cs_inv, pt);
      });
    }

  class distance_field
    {
    public:
      distance_field();
      distance_field(
        uint32_t i_x,
        uint32_t i_y,
        uint32_t i_z,
        float i_offset_x,
        float i_offset_y,
        float i_offset_z,
        const std::vector<vec3<float>>& vertices,
        const std::vector<vec3<uint32_t>>& triangles,
        bool signed_distance);

      distance_field& operator = (const distance_field&) = delete;

      void init(
        uint32_t i_x,
        uint32_t i_y,
        uint32_t i_z,
        float i_offset_x,
        float i_offset_y,
        float i_offset_z,
        const std::vector<vec3<float>>& vertices,
        const std::vector<vec3<uint32_t>>& triangles,
        bool signed_distance);
      void refresh_distance_field();
      void reset_offset(float i_offset_x, float i_offset_y, float i_offset_z);
      const uint32_t size_x() const;
      const uint32_t size_y() const;
      const uint32_t size_z() const;
      float offset_x() const;
      float offset_y() const;
      float offset_z() const;
      float distance_from_point_to_surface(const vec3<float>& i_point) const;
      boundingbox3d<float> bounding_box() const;

    private:
      uint32_t m_x;
      uint32_t m_y;
      uint32_t m_z;
      float m_offset_x;
      float m_offset_y;
      float m_offset_z;
      float m_step_x;
      float m_step_y;
      float m_step_z;
      std::vector<vec3<float>> m_vertices;
      std::vector<vec3<uint32_t>> m_triangles;
      std::vector<vec3<float>> triangle_normals;
      std::unique_ptr<qbvh> mp_qbvh;
      boundingbox3d<float> m_bounding_box;
      bool m_signed_distance;
    };

  class warping_tool
    {
    public:
      warping_tool(const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, int distance_field_discretization, float decay_factor, bool signed_distance);
      ~warping_tool();

      void transform(vec3<float>& pt);

      bool influences_point(const vec3<float>& pt) const;

      void set_decay_factor(float decay_factor);

      void set_rotation(const vec3<float>& center, const vec3<float>& vector, float angle);
      void set_translation(const vec3<float>& translation);

      const float4x4& get_coordinate_system() const;
      void set_coordinate_system(const float4x4& cs);

      void set_weight_type(int w);

    private:
      void _compute_boundingbox();
      void _transform(vec3<float>& pt, float distance) const;

    private:
      std::shared_ptr<distance_field> df;
      float4x4 coordinate_system, coordinate_system_inv;
      float4x4 Q, P, I;
      float decay_factor;
      boundingbox3d<float> bb;
      float rotation_angle;
      vec3<float> rotation_center;
      vec3<float> translation;
      int weight_type;
      bool signed_distance;
    };

  class pushpull_tool
    {
    public:
      pushpull_tool(const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, float decay_factor, bool signed_distance);
      ~pushpull_tool();

      void transform(vec3<float>& pt);

      bool influences_point(const vec3<float>& pt) const;

      void set_decay_factor(float decay_factor);

      void set_translation(const vec3<float>& translation);

      const float4x4& get_coordinate_system() const;
      void set_coordinate_system(const float4x4& cs);

      void set_weight_type(int w);

    private:
      void _compute_boundingbox();
      void _transform(vec3<float>& pt, float distance) const;

    private:
      std::unique_ptr<qbvh> mp_qbvh;
      std::vector<vec3<float>> m_vertices;
      std::vector<vec3<uint32_t>> m_triangles;
      std::vector<vec3<float>> triangle_normals;
      float4x4 coordinate_system, coordinate_system_inv;
      float decay_factor;
      boundingbox3d<float> bb;
      vec3<float> translation, transformed_translation;
      int weight_type;
      bool signed_distance;
    };
  }

#endif // JTK_DEFORMATION_H

#ifdef JTK_DEFORMATION_IMPLEMENTATION

namespace jtk
  {
  distance_field::distance_field()
    {
    }

  distance_field::distance_field(uint32_t i_x, uint32_t i_y, uint32_t i_z, float i_offset_x, float i_offset_y, float i_offset_z,
    const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, bool signed_distance)
    {
    init(i_x, i_y, i_z, i_offset_x, i_offset_y, i_offset_z, vertices, triangles, signed_distance);
    }

  void distance_field::init(uint32_t i_x, uint32_t i_y, uint32_t i_z, float i_offset_x, float i_offset_y, float i_offset_z,
    const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, bool signed_distance)
    {
    m_x = i_x;
    m_y = i_y;
    m_z = i_z;
    m_offset_x = i_offset_x;
    m_offset_y = i_offset_y;
    m_offset_z = i_offset_z;
    m_vertices = vertices;
    m_triangles = triangles;
    m_signed_distance = signed_distance;
    assert(m_x > 2);
    assert(m_y > 2);
    assert(m_z > 2);
    refresh_distance_field();
    }

  void distance_field::refresh_distance_field()
    {
    if (m_vertices.empty() || m_triangles.empty())
      return;
    compute_triangle_normals(triangle_normals, m_vertices.data(), m_triangles.data(), (uint32_t)m_triangles.size());
    m_bounding_box = bounding_volume_3d<float>(m_vertices.begin(), m_vertices.end());
    m_bounding_box.min[0] -= m_offset_x;
    m_bounding_box.min[1] -= m_offset_y;
    m_bounding_box.min[2] -= m_offset_z;
    m_bounding_box.max[0] += m_offset_x;
    m_bounding_box.max[1] += m_offset_y;
    m_bounding_box.max[2] += m_offset_z;

    m_step_x = (m_bounding_box.max[0] - m_bounding_box.min[0]) / (static_cast<float>(m_x) - float(1.0));
    m_step_y = (m_bounding_box.max[1] - m_bounding_box.min[1]) / (static_cast<float>(m_y) - float(1.0));
    m_step_z = (m_bounding_box.max[2] - m_bounding_box.min[2]) / (static_cast<float>(m_z) - float(1.0));

    mp_qbvh.reset(new qbvh(m_triangles, m_vertices.data()));
    }

  void distance_field::reset_offset(float i_offset_x, float i_offset_y, float i_offset_z)
    {
    m_offset_x = i_offset_x;
    m_offset_y = i_offset_y;
    m_offset_z = i_offset_z;
    refresh_distance_field();
    }

  const uint32_t distance_field::size_x() const
    {
    return m_x;
    }

  const uint32_t distance_field::size_y() const
    {
    return m_y;
    }

  const uint32_t distance_field::size_z() const
    {
    return m_z;
    }

  float distance_field::offset_x() const
    {
    return m_offset_x;
    }

  float distance_field::offset_y() const
    {
    return m_offset_y;
    }

  float distance_field::offset_z() const
    {
    return m_offset_z;
    }

  float distance_field::distance_from_point_to_surface(const vec3<float>& i_point) const
    {
    uint32_t triangle_id;
    auto h = mp_qbvh->find_closest_triangle(triangle_id, i_point, m_triangles.data(), m_vertices.data());
    if (m_signed_distance)
      {
      auto n = triangle_normals[triangle_id];
      const vec3<float> V0 = (m_vertices)[(m_triangles)[triangle_id][0]];
      const vec3<float> V1 = (m_vertices)[(m_triangles)[triangle_id][1]];
      const vec3<float> V2 = (m_vertices)[(m_triangles)[triangle_id][2]];
      const auto pos = V0 * (1.f - h.u - h.v) + h.u * V1 + h.v * V2;
      auto v = (i_point - pos);
      float sgn = dot(v, n);
      float dist = h.distance;
      return sgn >= 0.f ? dist : -dist;
      }
    else
      return h.distance;
    }


  boundingbox3d<float> distance_field::bounding_box() const
    {
    return m_bounding_box;
    }


  namespace
    {
    float _weight_very_soft(float i_distance, float i_decay_factor)
      {
      if (i_distance >= i_decay_factor)
        return 0.f;
      if (i_distance <= 0.f)
        return 1.f;
      float f = i_distance / i_decay_factor;
      f *= f;
      f -= 1.f;
      return f * f;
      }

    float _weight_soft(float i_distance, float i_decay_factor)
      {
      const float threshold = 0.75f;
      if (i_distance >= i_decay_factor)
        return 0.f;
      if (i_distance <= 0.f)
        return 1.f;
      if (i_distance >= i_decay_factor * (1.f - threshold))
        {
        float f = (i_distance - i_decay_factor * (1.f - threshold)) / (threshold * i_decay_factor);
        f *= f;
        f -= 1.f;
        return f * f;
        }
      return 1.f;
      }

    float _weight_hard(float i_distance, float i_decay_factor)
      {
      const float threshold = 0.5f;
      if (i_distance >= i_decay_factor)
        return 0.f;
      if (i_distance <= 0.f)
        return 1.f;
      if (i_distance >= i_decay_factor * (1.f - threshold))
        {
        float f = (i_distance - i_decay_factor * (1.f - threshold)) / (threshold * i_decay_factor);
        f *= f;
        f -= 1.f;
        return f * f;
        }
      return 1.f;
      }

    float _weight_very_hard(float i_distance, float i_decay_factor)
      {
      const float threshold = 0.25f;
      if (i_distance >= i_decay_factor)
        return 0.f;
      if (i_distance <= 0.f)
        return 1.f;
      if (i_distance >= i_decay_factor * (1.f - threshold))
        {
        float f = (i_distance - i_decay_factor * (1.f - threshold)) / (threshold * i_decay_factor);
        f *= f;
        f -= 1.f;
        return f * f;
        }
      return 1.f;
      }

    float _weight(float i_distance, float i_decay_factor, int type)
      {
      switch (type)
        {
        case 0: return _weight_very_soft(i_distance, i_decay_factor);
        case 1: return _weight_soft(i_distance, i_decay_factor);
        case 2: return _weight_hard(i_distance, i_decay_factor);
        case 3: return _weight_very_hard(i_distance, i_decay_factor);
        default: return _weight_very_soft(i_distance, i_decay_factor);
        }
      }
    }

  warping_tool::warping_tool(const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, int distance_field_discretization, float i_decay_factor, bool i_signed_distance) :
    rotation_angle(0.0), decay_factor(i_decay_factor), weight_type(0), signed_distance(i_signed_distance)
    {
    df = std::make_shared<distance_field>(distance_field_discretization, distance_field_discretization, distance_field_discretization, decay_factor * 1.2f, decay_factor * 1.2f, decay_factor * 1.2f, vertices, triangles, signed_distance);
    rotation_center = vec3<float>(0, 0, 0);
    translation = vec3<float>(0, 0, 0);
    coordinate_system = get_identity();
    coordinate_system_inv = get_identity();
    I = get_identity();
    _compute_boundingbox();
    }

  warping_tool::~warping_tool()
    {
    }

  void warping_tool::set_weight_type(int w)
    {
    weight_type = w;
    }

  void warping_tool::_transform(vec3<float>& pt, float distance) const
    {
    float w = _weight(distance, decay_factor, weight_type);
    if (w)
      {
      if (rotation_angle)
        {
        auto rot = P + (I - P) * std::cos(w * rotation_angle) + Q * std::sin(w * rotation_angle);
        pt = jtk::transform(rot, pt - rotation_center) + rotation_center + translation * w;
        }
      else
        pt = pt + translation * w;
      }
    }

  void warping_tool::transform(vec3<float>& pt)
    {
    vec3<float> transformed_pt = jtk::transform(coordinate_system_inv, pt);
    float distance = df->distance_from_point_to_surface(transformed_pt);
    _transform(pt, distance);
    }

  bool warping_tool::influences_point(const vec3<float>& pt) const
    {
    return inside(bb, pt);
    }

  void warping_tool::set_decay_factor(float decay)
    {
    decay_factor = decay;
    if (decay * 1.2 > df->offset_x())
      df->reset_offset(decay * 1.2f, decay * 1.2f, decay * 1.2f);
    _compute_boundingbox();
    }

  void warping_tool::set_rotation(const vec3<float>& center, const vec3<float>& vector, float angle)
    {
    rotation_angle = angle;
    rotation_center = center;
    Q = get_identity();
    P = get_identity();
    Q[0 + 4 * 0] = 0.f;
    Q[1 + 4 * 0] = vector[2];
    Q[2 + 4 * 0] = -vector[1];
    Q[0 + 4 * 1] = -vector[2];
    Q[1 + 4 * 1] = 0.f;
    Q[2 + 4 * 1] = vector[0];
    Q[0 + 4 * 2] = vector[1];
    Q[1 + 4 * 2] = -vector[0];
    Q[2 + 4 * 2] = 0.f;
    P[0 + 4 * 0] = vector[0] * vector[0];
    P[0 + 4 * 1] = vector[1] * vector[0];
    P[0 + 4 * 2] = vector[2] * vector[0];
    P[1 + 4 * 0] = vector[0] * vector[1];
    P[1 + 4 * 1] = vector[1] * vector[1];
    P[1 + 4 * 2] = vector[2] * vector[1];
    P[2 + 4 * 0] = vector[0] * vector[2];
    P[2 + 4 * 1] = vector[1] * vector[2];
    P[2 + 4 * 2] = vector[2] * vector[2];
    }

  void warping_tool::set_translation(const vec3<float>& t)
    {
    translation = t;
    }

  const float4x4& warping_tool::get_coordinate_system() const
    {
    return coordinate_system;
    }

  void warping_tool::set_coordinate_system(const float4x4& cs)
    {
    coordinate_system = cs;
    coordinate_system_inv = invert_orthonormal(cs);
    _compute_boundingbox();
    }

  void warping_tool::_compute_boundingbox()
    {
    bb = df->bounding_box();
    vec3<float> p[8];

    p[0] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.min[1], bb.min[2]));
    p[1] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.min[1], bb.min[2]));
    p[2] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.max[1], bb.min[2]));
    p[3] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.max[1], bb.min[2]));
    p[4] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.min[1], bb.max[2]));
    p[5] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.min[1], bb.max[2]));
    p[6] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.max[1], bb.max[2]));
    p[7] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.max[1], bb.max[2]));

    bb.min = p[0];
    bb.max = p[0];
    for (size_t i = 1; i < 8; ++i)
      add_point(bb, p[i]);
    }


  pushpull_tool::pushpull_tool(const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles, float i_decay_factor, bool i_signed_distance) :
    m_vertices(vertices), m_triangles(triangles), decay_factor(i_decay_factor), weight_type(0), signed_distance(i_signed_distance)
    {
    compute_triangle_normals(triangle_normals, m_vertices.data(), m_triangles.data(), (uint32_t)m_triangles.size());
    mp_qbvh.reset(new qbvh(m_triangles, m_vertices.data()));
    translation = vec3<float>(0, 1, 0);
    transformed_translation = translation;
    coordinate_system = get_identity();
    coordinate_system_inv = get_identity();
    _compute_boundingbox();
    }

  pushpull_tool::~pushpull_tool()
    {
    }

  void pushpull_tool::set_weight_type(int w)
    {
    weight_type = w;
    }

  void pushpull_tool::_transform(vec3<float>& pt, float distance) const
    {
    float w = _weight(distance, decay_factor, weight_type);
    if (w)
      {
      pt = pt + translation * w;
      }
    }

  void pushpull_tool::transform(vec3<float>& pt)
    {
    vec3<float> transformed_pt = jtk::transform(coordinate_system_inv, pt);

    uint32_t triangle_id;
    ray r;
    r.dir[0] = -transformed_translation[0];
    r.dir[1] = -transformed_translation[1];
    r.dir[2] = -transformed_translation[2];
    r.dir[3] = 0.f;
    r.orig[0] = transformed_pt[0];
    r.orig[1] = transformed_pt[1];
    r.orig[2] = transformed_pt[2];
    r.orig[3] = 1.f;
    r.t_near = -decay_factor;
    r.t_far = decay_factor;
    auto h = mp_qbvh->find_closest_triangle(triangle_id, r, m_triangles.data(), m_vertices.data());
    if (h.found)
      {
      _transform(pt, signed_distance ? h.distance : std::abs(h.distance));
      }
    }

  bool pushpull_tool::influences_point(const vec3<float>& pt) const
    {
    return inside(bb, pt);
    }

  void pushpull_tool::set_decay_factor(float decay)
    {
    decay_factor = decay;
    _compute_boundingbox();
    }

  void pushpull_tool::set_translation(const vec3<float>& t)
    {
    translation = t;
    transformed_translation = normalize(jtk::transform(coordinate_system_inv, translation, true));
    }

  const float4x4& pushpull_tool::get_coordinate_system() const
    {
    return coordinate_system;
    }

  void pushpull_tool::set_coordinate_system(const float4x4& cs)
    {
    coordinate_system = cs;
    coordinate_system_inv = invert_orthonormal(cs);
    transformed_translation = normalize(jtk::transform(coordinate_system_inv, translation, true));
    _compute_boundingbox();
    }

  void pushpull_tool::_compute_boundingbox()
    {
    bb = bounding_volume_3d<float>(m_vertices.begin(), m_vertices.end());
    bb.min[0] -= decay_factor;
    bb.min[1] -= decay_factor;
    bb.min[2] -= decay_factor;
    bb.max[0] += decay_factor;
    bb.max[1] += decay_factor;
    bb.max[2] += decay_factor;

    vec3<float> p[8];

    p[0] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.min[1], bb.min[2]));
    p[1] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.min[1], bb.min[2]));
    p[2] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.max[1], bb.min[2]));
    p[3] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.max[1], bb.min[2]));
    p[4] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.min[1], bb.max[2]));
    p[5] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.min[1], bb.max[2]));
    p[6] = jtk::transform(coordinate_system, vec3<float>(bb.min[0], bb.max[1], bb.max[2]));
    p[7] = jtk::transform(coordinate_system, vec3<float>(bb.max[0], bb.max[1], bb.max[2]));

    bb.min = p[0];
    bb.max = p[0];
    for (size_t i = 1; i < 8; ++i)
      add_point(bb, p[i]);
    }
  } // namespace jtk
#endif //JTK_DEFORMATION_IMPLEMENTATION