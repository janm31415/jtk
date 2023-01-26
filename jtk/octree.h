#ifndef JTK_OCTREE_H
#define JTK_OCTREE_H

#include "containers.h"
#include "vec.h"

namespace jtk
  {

  // single threaded octree
  template <class T>
  struct octree_traits
    {
    static const int32_t internal_bytes = sizeof(uint8_t);
    static const int32_t leaf_bytes = sizeof(T);
    static const int32_t pointer_bytes = sizeof(void*);

    typedef memory_pool<internal_bytes + 0 * pointer_bytes> _internal_node_with_0_children_pool;
    typedef memory_pool<internal_bytes + 1 * pointer_bytes> _internal_node_with_1_children_pool;
    typedef memory_pool<internal_bytes + 2 * pointer_bytes> _internal_node_with_2_children_pool;
    typedef memory_pool<internal_bytes + 3 * pointer_bytes> _internal_node_with_3_children_pool;
    typedef memory_pool<internal_bytes + 4 * pointer_bytes> _internal_node_with_4_children_pool;
    typedef memory_pool<internal_bytes + 5 * pointer_bytes> _internal_node_with_5_children_pool;
    typedef memory_pool<internal_bytes + 6 * pointer_bytes> _internal_node_with_6_children_pool;
    typedef memory_pool<internal_bytes + 7 * pointer_bytes> _internal_node_with_7_children_pool;
    typedef memory_pool<internal_bytes + 8 * pointer_bytes> _internal_node_with_8_children_pool;
    typedef memory_pool<leaf_bytes> _leaf_nodes_pool;
    };


  // use these traits for a concurrent octree
  template <class T>
  struct concurrent_octree_traits
    {
    static const int32_t internal_bytes = sizeof(uint8_t);
    static const int32_t leaf_bytes = sizeof(T);
    static const int32_t pointer_bytes = sizeof(void*);

    typedef concurrent_memory_pool<internal_bytes + 0 * pointer_bytes> _internal_node_with_0_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 1 * pointer_bytes> _internal_node_with_1_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 2 * pointer_bytes> _internal_node_with_2_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 3 * pointer_bytes> _internal_node_with_3_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 4 * pointer_bytes> _internal_node_with_4_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 5 * pointer_bytes> _internal_node_with_5_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 6 * pointer_bytes> _internal_node_with_6_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 7 * pointer_bytes> _internal_node_with_7_children_pool;
    typedef concurrent_memory_pool<internal_bytes + 8 * pointer_bytes> _internal_node_with_8_children_pool;
    typedef concurrent_memory_pool<leaf_bytes> _leaf_nodes_pool;
    };

  template <class T, class Traits = octree_traits<T>>
  class indexed_octree
    {
    public:
      static const int32_t internal_bytes = Traits::internal_bytes;
      static const int32_t leaf_bytes = Traits::leaf_bytes;
      static const int32_t pointer_bytes = Traits::pointer_bytes;

      indexed_octree(uint32_t max_depth) : _max_depth(max_depth)
        {
        _build_table();
        _root = _create_internal(0);
        }

      size_t memory_used() const
        {
        const size_t mem_int0 = _internal_0.memory_used();
        const size_t mem_int1 = _internal_1.memory_used();
        const size_t mem_int2 = _internal_2.memory_used();
        const size_t mem_int3 = _internal_3.memory_used();
        const size_t mem_int4 = _internal_4.memory_used();
        const size_t mem_int5 = _internal_5.memory_used();
        const size_t mem_int6 = _internal_6.memory_used();
        const size_t mem_int7 = _internal_7.memory_used();
        const size_t mem_int8 = _internal_8.memory_used();
        const size_t mem_leaf = _leaf_node.memory_used();
        return mem_int0 + mem_int1 + mem_int2 + mem_int3 + mem_int4 + mem_int5 + mem_int6 + mem_int7 + mem_int8 + mem_leaf;
        }

      bool in_bounds(const uint32_t* coord, uint32_t depth) const
        {
        if ((coord[0] >> (depth - 1)) > 1)
          return false;
        if ((coord[1] >> (depth - 1)) > 1)
          return false;
        if ((coord[2] >> (depth - 1)) > 1)
          return false;
        return true;
        }

      uint8_t* get_node(const uint32_t* coord, uint32_t depth) const
        {
        assert(in_bounds(coord, depth));
        uint8_t* node = _root;
        uint32_t idx = 0;
        for (uint32_t i = 1; i <= depth; ++i)
          {
          const uint32_t x = (coord[0] >> (depth - i)) & 1;
          const uint32_t y = (coord[1] >> (depth - i)) & 1;
          const uint32_t z = (coord[2] >> (depth - i)) & 1;
          idx = (z << 2) | (y << 1) | x;
          node = get_child(node, get_child_count(node, idx));
          }
        return node;
        }

      uint8_t* get_child(const uint8_t* parent, int count) const
        {
        assert(count >= 0 && count < 8);
        const uint8_t* location = parent + internal_bytes + pointer_bytes * count;
        uint8_t* ptr = nullptr;
        memcpy(&ptr, location, pointer_bytes);
        return ptr;
        }

      bool has_child(const uint8_t* parent, int index) const
        {
        assert(index >= 0 && index < 8);
        return ((parent[0] >> index) & 1) != 0;
        }

      int32_t get_number_of_children(const uint8_t* parent) const
        {
        return _number_of_children_table[parent[0]];
        }

      int32_t get_child_count(const uint8_t* parent, int index) const
        {
        assert(index >= 0 && index < 8);
        return _children_count_table[parent[0]][index];
        }

      int32_t get_child_index(const uint8_t* parent, int count) const
        {
        assert(count >= 0 && count < 8);
        return _children_index_table[parent[0]][count];
        }

      T get_value(const uint8_t* leaf) const
        {
        T val;
        memcpy(&val, leaf, sizeof(T));
        return val;
        }

      void set_value(uint8_t* leaf, const T& value)
        {
        memcpy(leaf, &value, sizeof(T));
        }

      uint8_t* add_leaf(const uint32_t* coord)
        {
        assert(in_bounds(coord, _max_depth));
        return _make_leaf(coord);
        }

      void remove_leaf(const uint32_t* coord)
        {
        assert(in_bounds(coord, _max_depth));
        return _remove_leaf(coord);
        }

      uint8_t* get_leaf(const uint32_t* coord) const
        {
        assert(in_bounds(coord, _max_depth));
        uint8_t* node = _root;
        for (uint32_t i = 1; i <= _max_depth; ++i)
          {
          const uint32_t x = (coord[0] >> (_max_depth - i)) & 1;
          const uint32_t y = (coord[1] >> (_max_depth - i)) & 1;
          const uint32_t z = (coord[2] >> (_max_depth - i)) & 1;
          const uint32_t idx = (z << 2) | (y << 1) | x;
          node = get_child(node, get_child_count(node, idx));
          }
        return node;
        }

      uint8_t* find_leaf(uint32_t* coord) const
        {
        assert(in_bounds(coord, _max_depth));
        uint8_t* node = _root;
        for (uint32_t i = 1; i <= _max_depth; ++i)
          {
          const uint32_t x = (coord[0] >> (_max_depth - i)) & 1;
          const uint32_t y = (coord[1] >> (_max_depth - i)) & 1;
          const uint32_t z = (coord[2] >> (_max_depth - i)) & 1;
          const uint32_t idx = (z << 2) | (y << 1) | x;
          if (!has_child(node, idx))
            return nullptr;
          node = get_child(node, get_child_count(node, idx));
          }
        return node;
        }

      uint8_t* find_leaf_neighbour(const uint32_t* coord, const uint32_t* offset, uint8_t* leaf_parent) const
        {
        assert(in_bounds(coord, _max_depth));
        assert(offset[0] <= 1);
        assert(offset[1] <= 1);
        assert(offset[2] <= 1);
        const uint32_t x = (coord[0] & 1) + offset[0];
        const uint32_t y = (coord[1] & 1) + offset[1];
        const uint32_t z = (coord[2] & 1) + offset[2];
        if (x <= 1 && y <= 1 && z <= 1) // same parent
          {
          const uint32_t idx = (z << 2) | (y << 1) | x;
          if (!has_child(leaf_parent, idx))
            return nullptr;
          return get_child(leaf_parent, get_child_count(leaf_parent, idx));
          }
        else
          {
          uint32_t new_coord[3] = { coord[0] + offset[0], coord[1] + offset[1], coord[2] + offset[2] };
          return find_leaf(new_coord);
          }
        }

      uint8_t* find_parent(uint32_t depth, const uint32_t* coord, int32_t& count) const
        {
        assert(in_bounds(coord, depth));
        uint8_t* node = _root;
        uint8_t* previous = nullptr;
        uint32_t idx = 0;
        for (uint32_t i = 1; i <= depth; ++i)
          {
          const uint32_t x = (coord[0] >> (depth - i)) & 1;
          const uint32_t y = (coord[1] >> (depth - i)) & 1;
          const uint32_t z = (coord[2] >> (depth - i)) & 1;
          idx = (z << 2) | (y << 1) | x;
          previous = node;
          node = get_child(node, get_child_count(node, idx));
          }
        count = previous ? get_child_count(previous, idx) : 0;
        assert(count >= 0 && count < 8 && (previous == nullptr || count < get_number_of_children(previous)));
        return previous;
        }

      void for_each_leaf(const std::function< void(uint8_t* /*node*/, uint8_t* /*parent*/, uint32_t* /*coord*/) >& action) const
        {
        uint32_t coord[3] = { 0, 0, 0 };
        _for_each_leaf(_root, nullptr, 0, coord, action);
        }

      bool tree_is_consistent(uint8_t* node = nullptr, uint32_t* coord = nullptr, uint32_t depth = 0) const
        {
        if (depth == _max_depth)
          return true;
        uint32_t origin[3] = { 0, 0, 0 };
        if (!node)
          {
          node = _root;
          coord = &origin[0];
          }
        int32_t nr_of_children = get_number_of_children(node);
        for (int32_t i = 0; i < nr_of_children; ++i)
          {
          uint8_t* child = get_child(node, i);
          int32_t cnt;
          int32_t idx = get_child_index(node, i);
          const uint32_t x = idx & 1;
          const uint32_t y = (idx >> 1) & 1;
          const uint32_t z = (idx >> 2) & 1;
          uint32_t new_coords[3] = { 2 * coord[0] + x, 2 * coord[1] + y, 2 * coord[2] + z };
          uint8_t* my_parent = find_parent(depth + 1, new_coords, cnt);
          if (my_parent != node)
            return false;
          bool res = tree_is_consistent(child, new_coords, depth + 1);
          if (!res)
            return false;
          }
        return true;
        }

      uint8_t* get_root() const
        {
        return _root;
        }

      uint32_t max_depth() const
        {
        return _max_depth;
        }

    private:

      void _for_each_leaf(uint8_t* parent, uint8_t* parents_parent, uint32_t depth, uint32_t* coords, const std::function< void(uint8_t*, uint8_t*, uint32_t*) >& action) const
        {
        if (depth == _max_depth) // leaf
          {
          action(parent, parents_parent, coords);
          }
        else
          {
          int32_t nr_of_children = get_number_of_children(parent);
          for (int32_t i = 0; i < nr_of_children; ++i)
            {
            uint8_t* child = get_child(parent, i);
            int idx = get_child_index(parent, i);
            const uint32_t x = idx & 1;
            const uint32_t y = (idx >> 1) & 1;
            const uint32_t z = (idx >> 2) & 1;
            uint32_t new_coords[3] = { 2 * coords[0] + x, 2 * coords[1] + y, 2 * coords[2] + z };
            _for_each_leaf(child, parent, depth + 1, new_coords, action);
            }
          }
        }

      bool _find_leaf_with_parent(uint8_t*& leaf, uint8_t*& parent, uint32_t* coord) const
        {
        leaf = _root;
        for (uint32_t i = 1; i <= _max_depth; ++i)
          {
          const uint32_t x = (coord[0] >> (_max_depth - i)) & 1;
          const uint32_t y = (coord[1] >> (_max_depth - i)) & 1;
          const uint32_t z = (coord[2] >> (_max_depth - i)) & 1;
          const uint32_t idx = (z << 2) | (y << 1) | x;
          if (!has_child(leaf, idx))
            return false;
          parent = leaf;
          leaf = get_child(parent, get_child_count(parent, idx));
          }
        return true;
        }

      void _set_child(uint8_t* parent, int32_t count, uint8_t* child) const
        {
        assert(count >= 0 && count < 8);
        uint8_t* location = parent + internal_bytes + pointer_bytes * count;
        memcpy(location, &child, pointer_bytes);
        }

      //returns the parent node of the child, so the "new" version of input pointer node
      uint8_t* _add_child(uint8_t* node, uint32_t* coord, uint32_t depth, int32_t child_index, uint8_t* parent, int32_t child_count)
        {
        assert(child_index >= 0 && child_index < 8);
        assert(!has_child(node, child_index));
        uint8_t child_mask = node[0];
        uint8_t* children[8] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
        int32_t nr_of_old_children = get_number_of_children(node);
        for (int32_t i = 0; i < nr_of_old_children; ++i)
          {
          int32_t idx = _children_index_table[child_mask][i];
          memcpy(&(children[idx]), node + internal_bytes + pointer_bytes * i, pointer_bytes);
          }
        auto lam = [&] {
          int32_t cnt;
          uint8_t* my_parent = find_parent(depth, coord, cnt);
          return cnt == child_count && my_parent == parent;
          };
        assert(lam());
        assert((parent == nullptr && node == _root) || (get_child(parent, child_count) == node));
        uint8_t new_child_mask = child_mask | (1 << child_index);
        assert(new_child_mask != child_mask);
        uint8_t* child = nullptr;
        if (depth == _max_depth - 1)
          child = _create_leaf();
        else
          child = _create_internal(0);
        assert(children[child_index] == nullptr);
        children[child_index] = child;
        int nr_of_children = _number_of_children_table[new_child_mask];
        uint8_t* new_node = _create_internal(nr_of_children);
        new_node[0] = new_child_mask;
        int32_t j = 0;
        for (int32_t i = 0; i < 8; ++i)
          {
          if (children[i])
            {
            memcpy(new_node + internal_bytes + pointer_bytes * j, &(children[i]), pointer_bytes);
            ++j;
            }
          }
        _remove_internal(node, nr_of_old_children);
        if (parent)
          memcpy(parent + internal_bytes + pointer_bytes * child_count, &new_node, pointer_bytes);
        else
          _root = new_node;
        return new_node;
        }

      //returns the parent node of the child, so the "new" version of input pointer node
      uint8_t* _remove_child(uint8_t* node, const uint32_t* coord, uint32_t depth, int32_t child_index, uint8_t* parent, int32_t child_count)
        {
        assert(child_index >= 0 && child_index < 8);
        assert(has_child(node, child_index));
        uint8_t child_mask = node[0];
        uint8_t* children[8] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
        int32_t nr_of_old_children = get_number_of_children(node);
        for (int32_t i = 0; i < nr_of_old_children; ++i)
          {
          int32_t idx = _children_index_table[child_mask][i];
          memcpy(&(children[idx]), node + internal_bytes + pointer_bytes * i, pointer_bytes);
          }
        auto lam = [&] {
          int32_t cnt;
          uint8_t* my_parent = find_parent(depth, coord, cnt);
          return cnt == child_count && my_parent == parent;
          };
        assert(lam());
        assert((parent == nullptr && node == _root) || (get_child(parent, child_count) == node));
        uint8_t new_child_mask = child_mask & ~(1 << child_index);
        assert(new_child_mask != child_mask);
        uint8_t* child = children[child_index];
        assert(child != nullptr);
        if (depth == _max_depth - 1)
          {
          _remove_leaf(child);
          }
        else
          {
          assert(get_number_of_children(child) == 0);
          _remove_internal(child, get_number_of_children(child));
          }
        children[child_index] = nullptr;
        int nr_of_children = _number_of_children_table[new_child_mask];
        uint8_t* new_node = _create_internal(nr_of_children);
        new_node[0] = new_child_mask;
        int32_t j = 0;
        for (int32_t i = 0; i < 8; ++i)
          {
          if (children[i])
            {
            memcpy(new_node + internal_bytes + pointer_bytes * j, &(children[i]), pointer_bytes);
            ++j;
            }
          }
        _remove_internal(node, nr_of_old_children);
        if (parent)
          memcpy(parent + internal_bytes + pointer_bytes * child_count, &new_node, pointer_bytes);
        else
          _root = new_node;
        return new_node;
        }

      uint8_t* _make_leaf(const uint32_t* coord)
        {
        uint8_t* node = _root;
        uint8_t* parent = nullptr;
        int32_t child_count = 0;
        for (uint32_t i = 1; i <= _max_depth; ++i)
          {
          const uint32_t x = (coord[0] >> (_max_depth - i)) & 1;
          const uint32_t y = (coord[1] >> (_max_depth - i)) & 1;
          const uint32_t z = (coord[2] >> (_max_depth - i)) & 1;
          const uint32_t idx = (z << 2) | (y << 1) | x;
          if (!has_child(node, idx))
            {
            uint32_t current_coord[3] = { coord[0] >> (_max_depth - i + 1), coord[1] >> (_max_depth - i + 1), coord[2] >> (_max_depth - i + 1) };
            node = _add_child(node, current_coord, i - 1, idx, parent, child_count);
            }
          assert(has_child(node, idx));
          parent = node;
          child_count = get_child_count(node, idx);
          node = get_child(node, child_count);
          }
        return node;
        }

      void _remove_leaf(const uint32_t* coord)
        {
        uint8_t* node = _root;
        uint8_t* parent = nullptr;
        int32_t child_count = 0;
        for (uint32_t i = 1; i < _max_depth; ++i)
          {
          const uint32_t x = (coord[0] >> (_max_depth - i)) & 1;
          const uint32_t y = (coord[1] >> (_max_depth - i)) & 1;
          const uint32_t z = (coord[2] >> (_max_depth - i)) & 1;
          const uint32_t idx = (z << 2) | (y << 1) | x;
          if (!has_child(node, idx))
            return; // nothing to remove          
          parent = node;
          child_count = get_child_count(node, idx);
          node = get_child(node, child_count);
          }
        const uint32_t x = coord[0] & 1;
        const uint32_t y = coord[1] & 1;
        const uint32_t z = coord[2] & 1;
        const uint32_t idx = (z << 2) | (y << 1) | x;
        if (!has_child(node, idx))
          return; // nothing to remove  
        uint32_t current_coord[3] = { coord[0] >> 1, coord[1] >> 1, coord[2] >> 1 };
        node = _remove_child(node, current_coord, _max_depth - 1, idx, parent, child_count);
        uint32_t depth = _max_depth - 1;
        while (parent != nullptr && get_number_of_children(node) == 0) // remove parent nodes as long as they are empty
          {
          uint32_t cidx = get_child_index(parent, child_count);
          node = parent;
          uint32_t cc[3] = { coord[0] >> (_max_depth - depth + 1), coord[1] >> (_max_depth - depth + 1), coord[2] >> (_max_depth - depth + 1) };
          parent = find_parent(depth - 1, cc, child_count);
          assert((parent == nullptr && node == _root) || get_child(parent, child_count) == node);
          node = _remove_child(node, cc, depth - 1, cidx, parent, child_count);
          --depth;
          }
        }

      inline uint8_t* _create_leaf()
        {
        uint8_t* node = _leaf_node.allocate();
        memset(node, 0, sizeof(T));
        return node;
        }

      inline uint8_t* _create_internal(int32_t number_of_children)
        {
        uint8_t* node = nullptr;
        switch (number_of_children)
          {
          case 0: node = _internal_0.allocate(); break;
          case 1: node = _internal_1.allocate(); break;
          case 2: node = _internal_2.allocate(); break;
          case 3: node = _internal_3.allocate(); break;
          case 4: node = _internal_4.allocate(); break;
          case 5: node = _internal_5.allocate(); break;
          case 6: node = _internal_6.allocate(); break;
          case 7: node = _internal_7.allocate(); break;
          case 8: node = _internal_8.allocate(); break;
          default: assert(0); break;
          }
        node[0] = 0;
        return node;
        }

      inline void _remove_leaf(uint8_t* node)
        {
        _leaf_node.deallocate(node);
        }

      inline void _remove_internal(uint8_t* node, int32_t number_of_children)
        {
        switch (number_of_children)
          {
          case 0: _internal_0.deallocate(node); break;
          case 1: _internal_1.deallocate(node); break;
          case 2: _internal_2.deallocate(node); break;
          case 3: _internal_3.deallocate(node); break;
          case 4: _internal_4.deallocate(node); break;
          case 5: _internal_5.deallocate(node); break;
          case 6: _internal_6.deallocate(node); break;
          case 7: _internal_7.deallocate(node); break;
          case 8: _internal_8.deallocate(node); break;
          default: assert(0); break;
          }
        }

      void _build_table()
        {
        for (int32_t i = 0; i < 256; ++i)
          {
          _number_of_children_table[i] = 0;
          int32_t count = 0;
          for (int32_t j = 0; j < 8; ++j)
            {
            _number_of_children_table[i] += ((i >> j) & 1);
            _children_count_table[i][j] = count;
            _children_index_table[i][count] = j;
            count += ((i >> j) & 1);
            }
          }
        }
    private:
      typename Traits::_internal_node_with_0_children_pool _internal_0;
      typename Traits::_internal_node_with_1_children_pool _internal_1;
      typename Traits::_internal_node_with_2_children_pool _internal_2;
      typename Traits::_internal_node_with_3_children_pool _internal_3;
      typename Traits::_internal_node_with_4_children_pool _internal_4;
      typename Traits::_internal_node_with_5_children_pool _internal_5;
      typename Traits::_internal_node_with_6_children_pool _internal_6;
      typename Traits::_internal_node_with_7_children_pool _internal_7;
      typename Traits::_internal_node_with_8_children_pool _internal_8;

      typename Traits::_leaf_nodes_pool _leaf_node;

      uint8_t* _root;
      uint32_t _max_depth;

      int32_t _number_of_children_table[256];
      int32_t _children_count_table[256][8];
      int32_t _children_index_table[256][8];
    };


  /*
  3d octree for points in the unit cube
  */
  template <class T, class Traits = octree_traits<T>>
  class octree
    {
    public:
      octree(uint32_t max_depth) : _oct(max_depth)
        {
        _dim = 1 << max_depth;
        }

      size_t memory_used() const
        {
        return _oct.memory_used();
        }

      template <class T2>
      uint8_t* add_leaf(const jtk::vec3<T2>& pt)
        {
        uint32_t coord[3];
        get_coordinate(coord, pt);
        return _oct.add_leaf(coord);
        }

      template <class T2>
      void remove_leaf(const jtk::vec3<T2>& pt)
        {
        uint32_t coord[3];
        get_coordinate(coord, pt);
        _oct.remove_leaf(coord);
        }

      template <class T2>
      uint8_t* get_leaf(const jtk::vec3<T2>& pt) const
        {
        uint32_t coord[3];
        get_coordinate(coord, pt);
        return _oct.get_leaf(coord);
        }

      template <class T2>
      uint8_t* find_leaf(const jtk::vec3<T2>& pt) const
        {
        uint32_t coord[3];
        get_coordinate(coord, pt);
        return _oct.find_leaf(coord);
        }

      bool tree_is_consistent() const
        {
        return _oct.tree_is_consistent();
        }

      uint32_t max_depth() const
        {
        return _oct.max_depth();
        }

      T get_value(const uint8_t* leaf) const
        {
        T val;
        memcpy(&val, leaf, sizeof(T));
        return val;
        }

      void set_value(uint8_t* leaf, const T& value)
        {
        memcpy(leaf, &value, sizeof(T));
        }

      void for_each_leaf(const std::function< void(uint8_t* /*leaf*/, const jtk::vec3<double>& /*leaf_center*/, double /*leaf_halfsize*/) >& action) const
        {
        const double leaf_halfsize = 0.5 / static_cast<double>(_dim);
        const double dim_inv = 1.0 / static_cast<double>(_dim);
        auto action2 = [&](uint8_t* node, uint8_t* /*parent*/, uint32_t* coord)
          {
          jtk::vec3<double> center(leaf_halfsize + static_cast<double>(coord[0]) * dim_inv, leaf_halfsize + static_cast<double>(coord[1]) * dim_inv, leaf_halfsize + static_cast<double>(coord[2]) * dim_inv);
          action(node, center, leaf_halfsize);
          };
        _oct.for_each_leaf(action2);
        }

    private:
      template <class T2>
      void get_coordinate(uint32_t* coord, const jtk::vec3<T2>& pt) const
        {
        coord[0] = static_cast<uint32_t>(std::trunc(pt[0] * static_cast<T2>(_dim)));
        coord[1] = static_cast<uint32_t>(std::trunc(pt[1] * static_cast<T2>(_dim)));
        coord[2] = static_cast<uint32_t>(std::trunc(pt[2] * static_cast<T2>(_dim)));
        }

    private:
      indexed_octree<T, Traits> _oct;
      uint32_t _dim;
    };
  }

#endif