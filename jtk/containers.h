#pragma once

#include <algorithm>
#include <list>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <stdexcept>

#include "concurrency.h"

namespace jtk
  {

  /////////////////////////////////////////
  // hashed_heap
  /////////////////////////////////////////

  namespace hashed_heap_details
    {

    template <typename _RandomAccessIterator, typename _Compare>
    inline void _up_heap(_RandomAccessIterator _first, _RandomAccessIterator _pos, _Compare _comp)
      {
      auto _index = _pos - _first;
      auto _parent = (_index - 1) / 2;
      auto _val = *_pos;

      while (_index > 0 && _comp(*(_first + _parent), _val)) {
        *(_first + _index) = *(_first + _parent);
        _index = _parent;
        _parent = (_parent - 1) / 2;
        }

      if (_pos != (_first + _index))
        *(_first + _index) = _val;
      }

    template <typename _RandomAccessIterator, typename _Compare>
    inline void _down_heap(_RandomAccessIterator _first,
      _RandomAccessIterator _last,
      _RandomAccessIterator _pos,
      _Compare _comp)
      {
      auto _len = _last - _first;
      auto _index = _pos - _first;
      auto _left = _index * 2 + 1;
      auto _right = _index * 2 + 2;
      auto _largest = _right;
      auto _val = *_pos;

      while (_index < _len) {
        if (_right < _len) {
          _largest = _comp(*(_first + _right), *(_first + _left)) ? _left : _right;
          }
        else if (_left < _len) {
          _largest = _left;
          }
        else {
          // Force termination
          _largest = _len;
          }

        if (_largest < _len && _comp(_val, *(_first + _largest))) {
          *(_first + _index) = *(_first + _largest);
          _index = _largest;
          _left = _index * 2 + 1;
          _right = _index * 2 + 2;
          }
        else
          break;
        }

      if (_pos != (_first + _index))
        *(_first + _index) = _val;
      }

    template <typename _RandomAccessIterator, typename _Compare>
    inline void _update_heap(_RandomAccessIterator _first,
      _RandomAccessIterator _last,
      _RandomAccessIterator _pos,
      _Compare _comp)
      {
      auto _index = (_pos - _first);
      auto _parent = (_index - 1) / 2;

      if (_index > 0 && _comp(*(_first + _parent), *(_pos)))
        _up_heap(_first, _pos, _comp);
      else
        _down_heap(_first, _last, _pos, _comp);
      }

    }

  template <class Key,
    class Data,
    class HashKey = std::hash<Key>,
    class KeyEquality = std::equal_to<Key>,
    class CompareData = std::less<Data>>
    class hashed_heap
    {
    private:
      using hash_type = std::unordered_map<Key, Data, HashKey, KeyEquality>;
      using heap_type = std::vector<Key>;

    public:
      using my_type = hashed_heap<Key, Data, HashKey, KeyEquality, CompareData>;
      using value_type = std::pair<Key, Data>;
      using const_reference = const std::pair<Key, Data>&;

    public:
      hashed_heap()
        : _hash()
        , _heap()
        {}

      hashed_heap(const hashed_heap& other)
        : _hash(other._hash)
        , _heap(other._heap)
        {}

      my_type& operator=(const my_type& right)
        {
        my_type temp(right);
        swap(temp);
        return (*this);
        }

      void swap(my_type& right)
        {
        _hash.swap(right._hash);
        _heap.swap(right._heap);
        }

      void reserve(std::size_t capacity) { _heap.reserve(capacity); }

      void push(const_reference val)
        {
        auto it = _hash.find(val.first);
        if (it != _hash.end() &&
          CompareData {}(val.second, it->second)) // element was already present, but with higher priority
          {
          it->second = val.second;
          auto iter =
            std::find_if(_heap.begin(), _heap.end(), [=](const Key& key) { return KeyEquality{}(key, val.first); });
          hashed_heap_details::_update_heap(_heap.begin(), _heap.end(), iter, Compare(_hash, _compare_data));
          }
        else if (it == _hash.end()) {
          _hash[val.first] = val.second;
          _heap.push_back(val.first);
          hashed_heap_details::_update_heap(_heap.begin(),
            _heap.end(),
            _heap.end() - 1,
            Compare(_hash, _compare_data)); // do not use std::push_heap, the STL does not guarantee that heaps
                                              // are implemented corresponding to _update_heap
          }
        }

      const hash_type& get_hash() const { return _hash; }

      const heap_type& get_heap() const { return _heap; }

      value_type top() const
        {
        std::pair<Key, Data> top_element(_heap.front(), _hash.at(_heap.front()));
        return top_element;
        }

      void pop() // do not use std::pop_heap, the STL does not guarantee that heaps are implemented corresponding to _update_heap
        {
        Key k = _heap.front();
        std::swap(_heap.front(), _heap.back());
        _heap.pop_back();
        if (!_heap.empty())
          hashed_heap_details::_update_heap(_heap.begin(), _heap.end(), _heap.begin(), Compare(_hash, _compare_data));
        _hash.erase(k);
        }

      bool empty() const { return _heap.empty(); }

    private:
      class Compare
        {
        private:
          hash_type& _hash;
          CompareData& _compare_data;

        public:

          Compare(const Compare& other)
            : _hash(other._hash)
            , _compare_data(other._compare_data)
            {}

          void swap(Compare& other)
            {
            std::swap(_hash, other._hash);
            std::swap(_compare_data, other._compare_data);
            }

          Compare& operator=(const Compare& other)
            {
            Compare temp(other);
            swap(temp);
            return *this;
            }

          Compare(hash_type& i_hash, CompareData& compare_data)
            : _hash(i_hash)
            , _compare_data(compare_data)
            {}

          bool operator()(const Key& left, const Key& right) const { return _compare_data(_hash[right], _hash[left]); }
        };

    private:
      hash_type _hash;
      heap_type _heap;
      CompareData _compare_data;
    };

  ///////////////////////////////////
  // memory_pool
  ///////////////////////////////////

  template <class T, int base = 16>
  class typed_memory_pool
    {
    private:
      size_t heap_size, heap_mask;
      std::vector<std::vector<T>> data;
      std::vector<std::vector<T*>> stack;
      size_t stack_size;
      size_t data_available;

    public:
      typed_memory_pool() : data(1), stack(1)
        {
        heap_size = 1 << base;
        heap_mask = heap_size - 1;

        data[0].resize(heap_size);
        stack[0].resize(heap_size);
        stack_size = heap_size;
        data_available = heap_size;

        for (uint32_t i = 0; i < heap_size; ++i)
          {
          stack[0][i] = data[0].data() + i;
          }
        }

      ~typed_memory_pool()
        {
        }

      void destroy()
        {
        data.clear();
        stack.clear();
        stack_size = 0;
        data_available = 0;
        }

      T* allocate()
        {
        if (data_available == 0)
          _allocate_data_block();
        --data_available;
        return stack[data_available >> base][data_available & heap_mask];
        }

      void deallocate(T* ptr)
        {
        if (data_available == stack_size)
          _allocate_stack_block();
        stack[data_available >> base][data_available & heap_mask] = ptr;
        ++data_available;
        }

      size_t allocated() const
        {
        return heap_size * data.size() - data_available;
        }

      size_t available() const
        {
        return data_available;
        }

      size_t size() const
        {
        return heap_size * data.size();
        }

      size_t memory_used() const
        {
        const size_t this_size = sizeof(typed_memory_pool<T, base>);
        const size_t outer_data_size = data.capacity() * sizeof(std::vector<T>);
        const size_t inner_data_size = data.size() * heap_size * sizeof(T);
        const size_t outer_stack_size = stack.capacity() * sizeof(std::vector<T*>);
        const size_t inner_stack_size = stack.size() * heap_size * sizeof(T*);
        return this_size + outer_data_size + inner_data_size + outer_stack_size + inner_stack_size;
        }

    private:
      void _allocate_data_block()
        {
        data.emplace_back(heap_size);
        for (uint32_t i = 0; i < heap_size; ++i)
          {
          stack[0][i] = data.back().data() + i;
          }
        data_available = heap_size;
        }

      void _allocate_stack_block()
        {
        stack.emplace_back(heap_size);
        stack_size += heap_size;
        }
    };

  template <class T, int base = 16>
  class concurrent_typed_memory_pool
    {
    private:
      size_t heap_size, heap_mask;
      std::vector<std::vector<T>> data;
      std::vector<std::vector<T*>> stack;
      size_t stack_size;
      std::atomic<size_t> data_available;
      spinlock_rw lock;

    public:
      concurrent_typed_memory_pool() : data(1), stack(1)
        {
        heap_size = 1 << base;
        heap_mask = heap_size - 1;

        data[0].resize(heap_size);
        stack[0].resize(heap_size);
        stack_size = heap_size;
        data_available = heap_size;

        for (uint32_t i = 0; i < heap_size; ++i)
          {
          stack[0][i] = data[0].data() + i;
          }
        }

      ~concurrent_typed_memory_pool()
        {
        }

      void destroy()
        {
        data.clear();
        stack.clear();
        stack_size = 0;
        data_available = 0;
        }

      T* allocate()
        {
        lock.lock_read();
        size_t expected = data_available.load();
        size_t desired = expected - 1;

        while (expected == 0)
          {
          lock.unlock();
          _allocate_data_block();
          lock.lock_read();
          expected = data_available.load();
          desired = expected - 1;
          }

        while (!data_available.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
          {
          expected = data_available.load();
          desired = expected - 1;

          while (expected == 0)
            {
            lock.unlock();
            _allocate_data_block();
            lock.lock_read();
            expected = data_available.load();
            desired = expected - 1;
            }
          }
        lock.unlock();
        return stack[desired >> base][desired & heap_mask];
        }

      void deallocate(T* ptr)
        {
        lock.lock_read();
        size_t expected = data_available.load();
        size_t desired = expected + 1;

        while (expected == stack_size)
          {
          lock.unlock();
          _allocate_stack_block();
          lock.lock_read();
          expected = data_available.load();
          desired = expected + 1;
          }

        while (!data_available.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
          {
          expected = data_available.load();
          desired = expected + 1;
          while (expected == stack_size)
            {
            lock.unlock();
            _allocate_stack_block();
            lock.lock_read();
            expected = data_available.load();
            desired = expected + 1;
            }
          }
        lock.unlock();
        stack[expected >> base][expected & heap_mask] = ptr;
        }

      size_t allocated() const
        {
        return heap_size * data.size() - data_available.load();
        }

      size_t available() const
        {
        return data_available.load();
        }

      size_t size() const
        {
        return heap_size * data.size();
        }

      size_t memory_used() const
        {
        const size_t this_size = sizeof(concurrent_typed_memory_pool<T, base>);
        const size_t outer_data_size = data.capacity() * sizeof(std::vector<T>);
        const size_t inner_data_size = data.size() * heap_size * sizeof(T);
        const size_t outer_stack_size = stack.capacity() * sizeof(std::vector<T*>);
        const size_t inner_stack_size = stack.size() * heap_size * sizeof(T*);
        return this_size + outer_data_size + inner_data_size + outer_stack_size + inner_stack_size;
        }

    private:
      void _allocate_data_block()
        {
        scoped_lock_rw scoped(lock, true);
        if (data_available.load() == 0)
          {
          data.emplace_back(heap_size);
          for (uint32_t i = 0; i < heap_size; ++i)
            {
            stack[0][i] = data.back().data() + i;
            }
          data_available = heap_size;
          }
        }

      void _allocate_stack_block()
        {
        scoped_lock_rw scoped(lock, true);
        if (data_available.load() == stack_size)
          {
          stack.emplace_back(heap_size);
          stack_size += heap_size;
          }
        }
    };


  template <int N, int base = 16>
  class memory_pool
    {
    private:
      size_t heap_size, heap_mask;
      uint8_t** data;
      uint8_t*** stack;
      size_t stack_size;
      size_t data_available;
      size_t number_of_data_blocks;
      size_t number_of_stack_blocks;

    public:
      memory_pool() : number_of_data_blocks(1), number_of_stack_blocks(1)
        {
        heap_size = 1 << base;
        heap_mask = heap_size - 1;

        data = (uint8_t**)malloc(sizeof(uint8_t*));
        data[0] = (uint8_t*)malloc(heap_size * N);
        stack = (uint8_t***)malloc(sizeof(uint8_t**));
        stack[0] = (uint8_t**)malloc(heap_size * sizeof(uint8_t*));
        stack_size = heap_size;
        data_available = heap_size;

        for (uint32_t i = 0; i < heap_size; ++i)
          {
          stack[0][i] = data[0] + i * N;
          }
        }

      ~memory_pool()
        {
        destroy();
        }

      void destroy()
        {
        for (size_t i = 0; i < number_of_data_blocks; ++i)
          {
          free(data[i]);
          }
        for (size_t i = 0; i < number_of_stack_blocks; ++i)
          {
          free(stack[i]);
          }
        free(data);
        free(stack);
        data = nullptr;
        stack = nullptr;
        stack_size = 0;
        data_available = 0;
        number_of_data_blocks = 0;
        number_of_stack_blocks = 0;
        }

      uint8_t* allocate()
        {
        if (data_available == 0)
          _allocate_data_block();
        --data_available;
        return stack[data_available >> base][data_available & heap_mask];
        }

      void deallocate(uint8_t* ptr)
        {
        if (data_available == stack_size)
          _allocate_stack_block();
        stack[data_available >> base][data_available & heap_mask] = ptr;
        ++data_available;
        }

      size_t allocated() const
        {
        return heap_size * number_of_data_blocks - data_available;
        }

      size_t available() const
        {
        return data_available;
        }

      size_t size() const
        {
        return heap_size * number_of_data_blocks;
        }

      size_t memory_used() const
        {
        const size_t this_size = sizeof(memory_pool<N, base>);
        const size_t inner_data_size = number_of_data_blocks * heap_size * N;
        const size_t inner_stack_size = number_of_stack_blocks * heap_size * sizeof(uint8_t*);
        return this_size + inner_data_size + inner_stack_size;
        }

    private:
      void _allocate_data_block()
        {
        ++number_of_data_blocks;
        data = (uint8_t**)realloc(data, sizeof(uint8_t*) * number_of_data_blocks);
        data[number_of_data_blocks - 1] = (uint8_t*)malloc(heap_size * N);
        for (uint32_t i = 0; i < heap_size; ++i)
          {
          stack[0][i] = data[number_of_data_blocks - 1] + i * N;
          }
        data_available = heap_size;
        }

      void _allocate_stack_block()
        {
        ++number_of_stack_blocks;
        stack_size += heap_size;
        stack = (uint8_t***)realloc(stack, sizeof(uint8_t**) * number_of_stack_blocks);
        stack[number_of_stack_blocks - 1] = (uint8_t**)malloc(heap_size * sizeof(uint8_t*));
        }
    };

  template <int N, int base = 16>
  class concurrent_memory_pool
    {
    private:
      size_t heap_size, heap_mask;
      uint8_t** data;
      uint8_t*** stack;
      size_t stack_size;
      std::atomic<size_t> data_available;
      size_t number_of_data_blocks;
      size_t number_of_stack_blocks;
      spinlock_rw lock;

    public:
      concurrent_memory_pool() : number_of_data_blocks(1), number_of_stack_blocks(1)
        {
        heap_size = 1 << base;
        heap_mask = heap_size - 1;

        data = (uint8_t**)malloc(sizeof(uint8_t*));
        data[0] = (uint8_t*)malloc(heap_size * N);
        stack = (uint8_t***)malloc(sizeof(uint8_t**));
        stack[0] = (uint8_t**)malloc(heap_size * sizeof(uint8_t*));
        stack_size = heap_size;
        data_available = heap_size;

        for (uint32_t i = 0; i < heap_size; ++i)
          {
          stack[0][i] = data[0] + i * N;
          }
        }

      ~concurrent_memory_pool()
        {
        destroy();
        }

      void destroy()
        {
        for (size_t i = 0; i < number_of_data_blocks; ++i)
          {
          free(data[i]);
          }
        for (size_t i = 0; i < number_of_stack_blocks; ++i)
          {
          free(stack[i]);
          }
        free(data);
        free(stack);
        data = nullptr;
        stack = nullptr;
        stack_size = 0;
        data_available = 0;
        number_of_data_blocks = 0;
        number_of_stack_blocks = 0;
        }

      uint8_t* allocate()
        {
        lock.lock_read();
        size_t expected = data_available.load();
        size_t desired = expected - 1;

        while (expected == 0)
          {
          lock.unlock();
          _allocate_data_block();
          lock.lock_read();
          expected = data_available.load();
          desired = expected - 1;
          }

        while (!data_available.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
          {
          expected = data_available.load();
          desired = expected - 1;

          while (expected == 0)
            {
            lock.unlock();
            _allocate_data_block();
            lock.lock_read();
            expected = data_available.load();
            desired = expected - 1;
            }
          }
        lock.unlock();
        return stack[desired >> base][desired & heap_mask];
        }

      void deallocate(uint8_t* ptr)
        {
        lock.lock_read();
        size_t expected = data_available.load();
        size_t desired = expected + 1;

        while (expected == stack_size)
          {
          lock.unlock();
          _allocate_stack_block();
          lock.lock_read();
          expected = data_available.load();
          desired = expected + 1;
          }

        while (!data_available.compare_exchange_weak(expected, desired, std::memory_order_release, std::memory_order_relaxed))
          {
          expected = data_available.load();
          desired = expected + 1;
          while (expected == stack_size)
            {
            lock.unlock();
            _allocate_stack_block();
            lock.lock_read();
            expected = data_available.load();
            desired = expected + 1;
            }
          }
        lock.unlock();
        stack[expected >> base][expected & heap_mask] = ptr;
        }

      size_t allocated() const
        {
        return heap_size * number_of_data_blocks - data_available.load();
        }

      size_t available() const
        {
        return data_available.load();
        }

      size_t size() const
        {
        return heap_size * number_of_data_blocks;
        }

      size_t memory_used() const
        {
        const size_t this_size = sizeof(concurrent_memory_pool<N, base>);
        const size_t inner_data_size = number_of_data_blocks * heap_size * N;
        const size_t inner_stack_size = number_of_stack_blocks * heap_size * sizeof(uint8_t*);
        return this_size + inner_data_size + inner_stack_size;
        }

    private:
      void _allocate_data_block()
        {
        scoped_lock_rw scoped(lock, true);
        if (data_available.load() == 0)
          {
          ++number_of_data_blocks;
          data = (uint8_t**)realloc(data, sizeof(uint8_t*) * number_of_data_blocks);
          data[number_of_data_blocks - 1] = (uint8_t*)malloc(heap_size * N);
          for (uint32_t i = 0; i < heap_size; ++i)
            {
            stack[0][i] = data[number_of_data_blocks - 1] + i * N;
            }
          data_available = heap_size;
          }
        }

      void _allocate_stack_block()
        {
        scoped_lock_rw scoped(lock, true);
        if (data_available.load() == stack_size)
          {
          ++number_of_stack_blocks;
          stack_size += heap_size;
          stack = (uint8_t***)realloc(stack, sizeof(uint8_t**) * number_of_stack_blocks);
          stack[number_of_stack_blocks - 1] = (uint8_t**)malloc(heap_size * sizeof(uint8_t*));
          }
        }
    };


  ///////////////////////////////////////
  // flat_map
  ///////////////////////////////////////

  template <typename Key, typename T, typename flat_map_type>
  class flat_map_iterator
    {
    public:
      typedef flat_map_iterator<Key, T, flat_map_type> self_type;
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef T value_type;
      typedef T& reference;
      typedef T* pointer;
      typedef std::ptrdiff_t difference_type;

      flat_map_iterator() : _si(nullptr)
        {
        }

      flat_map_iterator(flat_map_type* si, bool begin) : _si(si)
        {
        if (_si)
          {
          _it = begin ? _si->m.begin() : _si->m.end();
          }
        }

      reference operator*() const
        {
        return _it->second;
        }

      pointer operator -> () const
        {
        return &_it->second;
        }

      Key key() const
        {
        return _it->first;
        }

      self_type& operator++()
        {
        ++_it;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        --_it;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return (_it == other._it);
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:
      flat_map_type* _si;
      typename std::vector<typename flat_map_type::value_type>::iterator _it;
    };

  template <typename Key, typename T, typename flat_map_type>
  class flat_map_const_iterator
    {
    public:
      typedef flat_map_const_iterator<Key, T, flat_map_type> self_type;
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef T value_type;
      typedef const T& reference;
      typedef const T* pointer;
      typedef std::ptrdiff_t difference_type;


      flat_map_const_iterator() : _si(nullptr)
        {
        }

      flat_map_const_iterator(const flat_map_type* si, bool begin) : _si(si)
        {
        if (_si)
          {
          _it = begin ? _si->m.begin() : _si->m.end();
          }
        }

      reference operator*() const
        {
        return _it->second;
        }

      pointer operator -> () const
        {
        return &_it->second;
        }

      Key key() const
        {
        return _it->first;
        }

      self_type& operator++()
        {
        ++_it;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        --_it;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return (_it == other._it);
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:
      const flat_map_type* _si;
      typename std::vector<typename flat_map_type::value_type>::const_iterator _it;
    };

  template <typename Key, typename T>
  class flat_map
    {
    public:
      typedef flat_map<Key, T> self_type;
      typedef flat_map_iterator<Key, T, flat_map<Key, T>> iterator;
      typedef flat_map_const_iterator<Key, T, flat_map<Key, T>> const_iterator;
      typedef std::pair<Key, T> value_type;

      flat_map()
        {
        }

      ~flat_map()
        {
        }

      flat_map(const self_type& other) : m(other.m)
        {
        }

      flat_map(self_type&& other) : flat_map()
        {
        swap(other);
        }

      self_type& operator=(const self_type& other)
        {
        self_type temp(other);
        swap(temp);
        return *this;
        }

      self_type& operator = (self_type&& other)
        {
            {
            self_type empty_temp;
            swap(empty_temp);
            }
            swap(other);
            return *this;
        }

      void swap(self_type& other)
        {
        std::swap(m, other.m);
        }

      T& put(const Key& key)
        {
        return get_value(key).second;
        }

      const T& get(const Key& key) const
        {
        return get_value(key).second;
        }

      bool has(const Key& key) const
        {
        return has_key(key);
        }

      iterator begin()
        {
        return iterator(this, true);
        }

      iterator end()
        {
        return iterator(this, false);
        }

      const_iterator begin() const
        {
        return const_iterator(this, true);
        }

      const_iterator end() const
        {
        return const_iterator(this, false);
        }

      const_iterator cbegin() const
        {
        return const_iterator(this, true);
        }

      const_iterator cend() const
        {
        return const_iterator(this, false);
        }

      uint64_t memory_used() const
        {
        uint64_t mem = m.capacity() * sizeof(value_type);
        return mem;
        }

      uint64_t size() const
        {
        return m.size();
        }

      bool empty() const
        {
        return m.empty();
        }

      size_t capacity() const
        {
        return m.capacity();
        }

      void reserve(size_t sz)
        {
        m.reserve(sz);
        }

    private:
      inline value_type& get_value(const Key& key)
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        if (it == m.end() || it->first != key)
          {
          value_type new_entry;
          new_entry.first = key;
          it = m.insert(it, new_entry);
          }
        return *it;
        }

      inline const value_type& get_value(const Key& key) const
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        return *it;
        }

      inline bool has_key(const Key& key) const
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        return !(it == m.end() || it->first != key);
        }

    private:
      std::vector<value_type> m;

      template <typename OtherKey, typename OtherT, typename flat_map_type>
      friend class flat_map_const_iterator;

      template <typename OtherKey, typename OtherT, typename flat_map_type>
      friend class flat_map_iterator;
    };

  ////////////////////////////////
  // flat_list_map
  ////////////////////////////////

  template <typename Key, typename T, typename flat_list_map_type>
  class flat_list_map_iterator
    {
    public:
      typedef flat_list_map_iterator<Key, T, flat_list_map_type> self_type;
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef T value_type;
      typedef T& reference;
      typedef T* pointer;
      typedef std::ptrdiff_t difference_type;

      flat_list_map_iterator() : _si(nullptr)
        {
        }

      flat_list_map_iterator(flat_list_map_type* si, bool begin) : _si(si)
        {
        if (_si)
          {
          _it = begin ? _si->m.begin() : _si->m.end();
          }
        }

      reference operator*() const
        {
        return _it->second;
        }

      pointer operator -> () const
        {
        return &_it->second;
        }

      Key key() const
        {
        return _it->first;
        }

      self_type& operator++()
        {
        ++_it;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        --_it;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return (_it == other._it);
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:
      flat_list_map_type* _si;
      typename std::list<typename flat_list_map_type::value_type>::iterator _it;
    };

  template <typename Key, typename T, typename flat_list_map_type>
  class flat_list_map_const_iterator
    {
    public:
      typedef flat_list_map_const_iterator<Key, T, flat_list_map_type> self_type;
      typedef std::bidirectional_iterator_tag iterator_category;
      typedef T value_type;
      typedef const T& reference;
      typedef const T* pointer;
      typedef std::ptrdiff_t difference_type;


      flat_list_map_const_iterator() : _si(nullptr)
        {
        }

      flat_list_map_const_iterator(const flat_list_map_type* si, bool begin) : _si(si)
        {
        if (_si)
          {
          _it = begin ? _si->m.begin() : _si->m.end();
          }
        }

      reference operator*() const
        {
        return _it->second;
        }

      pointer operator -> () const
        {
        return &_it->second;
        }

      Key key() const
        {
        return _it->first;
        }

      self_type& operator++()
        {
        ++_it;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        --_it;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return (_it == other._it);
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:
      const flat_list_map_type* _si;
      typename std::list<typename flat_list_map_type::value_type>::const_iterator _it;
    };

  template <typename Key, typename T>
  class flat_list_map
    {
    public:
      typedef flat_list_map<Key, T> self_type;
      typedef flat_list_map_iterator<Key, T, flat_list_map<Key, T>> iterator;
      typedef flat_list_map_const_iterator<Key, T, flat_list_map<Key, T>> const_iterator;
      typedef std::pair<Key, T> value_type;

      flat_list_map()
        {
        }

      ~flat_list_map()
        {
        }

      flat_list_map(const self_type& other) : m(other.m)
        {
        }

      flat_list_map(self_type&& other) : flat_list_map()
        {
        swap(other);
        }

      self_type& operator=(const self_type& other)
        {
        self_type temp(other);
        swap(temp);
        return *this;
        }

      self_type& operator = (self_type&& other)
        {
            {
            self_type empty_temp;
            swap(empty_temp);
            }
            swap(other);
            return *this;
        }

      void swap(self_type& other)
        {
        std::swap(m, other.m);
        }

      T& put(const Key& key)
        {
        return get_value(key).second;
        }

      const T& get(const Key& key) const
        {
        return get_value(key).second;
        }

      bool has(const Key& key) const
        {
        return has_key(key);
        }

      void erase(const Key& key)
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        m.erase(it);
        }

      iterator begin()
        {
        return iterator(this, true);
        }

      iterator end()
        {
        return iterator(this, false);
        }

      const_iterator begin() const
        {
        return const_iterator(this, true);
        }

      const_iterator end() const
        {
        return const_iterator(this, false);
        }

      const_iterator cbegin() const
        {
        return const_iterator(this, true);
        }

      const_iterator cend() const
        {
        return const_iterator(this, false);
        }

      uint64_t memory_used() const
        {
        uint64_t mem = m.size() * sizeof(value_type);
        return mem;
        }

      uint64_t size() const
        {
        return m.size();
        }

      bool empty() const
        {
        return m.empty();
        }

      size_t capacity() const
        {
        return m.capacity();
        }

      void reserve(size_t sz)
        {
        m.reserve(sz);
        }

    private:
      inline value_type& get_value(const Key& key)
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        if (it == m.end() || it->first != key)
          {
          value_type new_entry;
          new_entry.first = key;
          it = m.insert(it, new_entry);
          }
        return *it;
        }

      inline const value_type& get_value(const Key& key) const
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        return *it;
        }

      inline bool has_key(const Key& key) const
        {
        value_type dummy(key, T());
        auto it = std::lower_bound(m.begin(), m.end(), dummy, [&](const value_type& left, const value_type& right)
          {
          return left.first < right.first;
          });
        return !(it == m.end() || it->first != key);
        }

    private:
      std::list<value_type> m;

      template <typename OtherKey, typename OtherT, typename flat_list_map_type>
      friend class flat_list_map_const_iterator;

      template <typename OtherKey, typename OtherT, typename flat_list_map_type>
      friend class flat_list_map_iterator;
    };


  ////////////////////////////////////////////
  // indexmap
  ////////////////////////////////////////////

  template <typename T, typename indexmap_type>
  class indexmap_const_iterator
    {
    public:
      typedef indexmap_const_iterator<T, indexmap_type> self_type;
      typedef std::forward_iterator_tag iterator_category;
      typedef T value_type;
      typedef const T& reference;
      typedef std::ptrdiff_t difference_type;
      typedef flat_map<uint64_t, T> bucket;

      indexmap_const_iterator() : _si(nullptr) {}

      indexmap_const_iterator(const indexmap_type* si) : _si(si)
        {
        if (_si)
          {
          bucket_it = _si->table.begin();
          entry_it = bucket_it->begin();
          if (entry_it == bucket_it->end())
            _find_next_nonempty_bucket();
          }
        }

      reference operator*() const
        {
        return *entry_it;
        }

      uint64_t entry() const
        {
        return entry_it.key();
        }

      self_type& operator++()
        {
        ++entry_it;
        if (entry_it == bucket_it->end())
          _find_next_nonempty_bucket();
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return ((_si == other._si) && (!_si || ((bucket_it == other.bucket_it) && (entry_it == other.entry_it))));
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:

      void _find_next_nonempty_bucket()
        {
        ++bucket_it;
        while (bucket_it != _si->table.end() && bucket_it->empty())
          ++bucket_it;
        if (bucket_it == _si->table.end())
          {
          _si = nullptr;
          }
        else
          {
          entry_it = bucket_it->begin();
          }
        }

    private:
      const indexmap_type* _si;
      typename std::vector<bucket>::const_iterator bucket_it;
      typename flat_map<uint64_t, T>::const_iterator entry_it;
    };

  template <class T>
  class indexmap
    {
    public:
      typedef indexmap<T> self_type;
      typedef indexmap_const_iterator<T, indexmap<T>> const_iterator;
      typedef flat_map<uint64_t, T> bucket;

      indexmap(uint32_t nr_of_buckets = 4096)
        {
        if (nr_of_buckets == 0)
          nr_of_buckets = 4096;
        nr_of_buckets = get_nearest_power_of_two(nr_of_buckets);
        table.resize(nr_of_buckets);
        mask = nr_of_buckets - 1;
        }

      ~indexmap()
        {
        }

      indexmap(const self_type& other) : table(other.table), mask(other.mask)
        {
        }

      indexmap(self_type&& other) : indexmap()
        {
        swap(other);
        }

      self_type& operator=(const self_type& other)
        {
        self_type temp(other);
        swap(temp);
        return *this;
        }

      self_type& operator = (self_type&& other)
        {
            {
            self_type empty_temp;
            swap(empty_temp);
            }
            swap(other);
            return *this;
        }

      void swap(self_type& other)
        {
        std::swap(table, other.table);
        std::swap(mask, other.mask);
        }

      T& put(uint64_t index)
        {
        bucket& b = get_bucket(index);
        return get_entry(index, b);
        }

      const T& get(uint64_t index) const
        {
        const bucket& b = get_bucket(index);
        return get_entry(index, b);
        }

      bool has(uint64_t index) const
        {
        const bucket& b = get_bucket(index);
        return has_entry(index, b);
        }

      const_iterator begin() const
        {
        return const_iterator(this);
        }

      const_iterator end() const
        {
        return const_iterator(nullptr);
        }

      const_iterator cbegin() const
        {
        return const_iterator(this);
        }

      const_iterator cend() const
        {
        return const_iterator(nullptr);
        }

      uint64_t memory_used() const
        {
        uint64_t mem = table.size() * sizeof(bucket);
        for (const auto& b : table)
          {
          mem += b.memory_used();
          }
        return mem;
        }

      uint64_t size() const
        {
        uint64_t _sz = 0;
        for (const auto& b : table)
          _sz += b.size();
        return _sz;
        }

      uint64_t bucket_count() const
        {
        return (uint64_t)table.size();
        }

      uint64_t bucket_size(uint32_t idx) const
        {
        return (uint64_t)table[idx].size();
        }

      double load_factor() const
        {
        return (double)size() / (double)bucket_count();
        }

    private:

      inline uint32_t get_nearest_power_of_two(uint32_t v)
        {
        --v;
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;
        ++v;
        return v;
        }

      inline T& get_entry(uint64_t index, bucket& b)
        {
        return b.put(index);
        }

      inline const T& get_entry(uint64_t index, const bucket& b) const
        {
        return b.get(index);
        }

      inline bool has_entry(uint64_t index, const bucket& b) const
        {
        return b.has(index);
        }

      inline const bucket& get_bucket(uint64_t index) const
        {
        return table[index & mask];
        }

      inline bucket& get_bucket(uint64_t index)
        {
        return const_cast<bucket&>(static_cast<const self_type&>(*this).get_bucket(index));
        }

    private:
      uint32_t mask;
      std::vector<bucket> table;

      template <typename OtherT, typename indexmap_type>
      friend class indexmap_const_iterator;
    };


  //////////////////////////////////////
  // sparse_buffer
  //////////////////////////////////////

  namespace sparse_buffer_details
    {

    template <int A>
    struct get_power_two
      {
      static const int value = 2 * get_power_two<A - 1>::value;
      };

    template <>
    struct get_power_two<0>
      {
      static const int value = 1;
      };
    }


  template <typename T, typename sparse_buffer_type>
  class sparse_buffer_const_iterator
    {
    public:
      typedef sparse_buffer_const_iterator<T, sparse_buffer_type> self_type;
      typedef std::forward_iterator_tag iterator_category;
      typedef T value_type;
      typedef const T& reference;
      typedef std::ptrdiff_t difference_type;
      typedef typename std::pair<uint64_t, T> bucket;

      sparse_buffer_const_iterator() : _si(nullptr) {}

      sparse_buffer_const_iterator(const sparse_buffer_type* si) : _si(si)
        {
        if (_si)
          {
          v_it = _si->v.begin();
          m_it = _si->m.begin();
          }
        }

      reference operator*() const
        {
        return (v_it == _si->v.end()) ? *m_it : v_it->second;
        }

      uint64_t entry() const
        {
        return (v_it == _si->v.end()) ? m_it.entry() : v_it->first;
        }

      self_type& operator++()
        {
        if (v_it == _si->v.end())
          {
          ++m_it;
          if (m_it == _si->m.end())
            _si = nullptr;
          }
        else
          {
          ++v_it;
          while (v_it != _si->v.end() && v_it->first == (uint64_t)-1)
            ++v_it;
          if (v_it == _si->v.end() && m_it == _si->m.end())
            _si = nullptr;
          }
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return ((_si == other._si) && (!_si || (v_it == other.v_it && m_it == other.m_it)));
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:
      const sparse_buffer_type* _si;
      typename std::vector<bucket>::const_iterator v_it;
      typename indexmap<T>::const_iterator m_it;
    };

  template <class T, int cluster_power = 1>
  class sparse_buffer
    {
    public:
      typedef sparse_buffer<T, cluster_power> self_type;
      typedef sparse_buffer_const_iterator<T, sparse_buffer<T, cluster_power>> const_iterator;
      typedef std::pair<uint64_t, T> value_type;

      sparse_buffer(uint32_t nr_of_buckets = 4096) : m(nr_of_buckets)
        {
        if (nr_of_buckets == 0)
          nr_of_buckets = 4096;
        nr_of_buckets = get_nearest_power_of_two(nr_of_buckets);
        value_type dummy((uint64_t)-1, T());
        v.resize(nr_of_buckets << cluster_power, dummy);
        mask = nr_of_buckets - 1;
        }

      ~sparse_buffer()
        {
        }

      sparse_buffer(const self_type& other) : v(other.v), mask(other.mask), m(other.m)
        {
        }

      sparse_buffer(self_type&& other) : sparse_buffer()
        {
        swap(other);
        }

      self_type& operator=(const self_type& other)
        {
        self_type temp(other);
        swap(temp);
        return *this;
        }

      self_type& operator = (self_type&& other)
        {
            {
            self_type empty_temp;
            swap(empty_temp);
            }
            swap(other);
            return *this;
        }

      void swap(self_type& other)
        {
        std::swap(v, other.v);
        m.swap(other.m);
        std::swap(mask, other.mask);
        }

      T& put(uint64_t index)
        {
        return get_entry(index);
        }

      const T& get(uint64_t index) const
        {
        return get_entry(index);
        }

      bool has(uint64_t index) const
        {
        return has_entry(index);
        }

      const_iterator begin() const
        {
        return const_iterator(this);
        }

      const_iterator end() const
        {
        return const_iterator(nullptr);
        }

      const_iterator cbegin() const
        {
        return const_iterator(this);
        }

      const_iterator cend() const
        {
        return const_iterator(nullptr);
        }

      uint64_t memory_used() const
        {
        return v.capacity() * sizeof(value_type) + m.memory_used();
        }

      uint64_t size() const
        {
        uint64_t _sz = 0;
        for (const auto& pr : v)
          {
          if (pr.first != (uint64_t)(-1))
            ++_sz;
          }
        _sz += m.size();
        return _sz;
        }

      uint64_t bucket_count() const
        {
        return (uint64_t)(mask + 1);
        }

      double load_factor() const
        {
        return (double)size() / (double)bucket_count();
        }

    private:

      inline uint32_t get_nearest_power_of_two(uint32_t n)
        {
        --n;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        ++n;
        return n;
        }

      inline T& get_entry(uint64_t index)
        {
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == (uint64_t)-1)
            {
            v[ind].first = index;
            return v[ind].second;
            }
          if (v[ind].first == index)
            return v[ind].second;
          }
        return m.put(index);
        }

      inline const T& get_entry(uint64_t index) const
        {
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == index)
            return v[ind].second;
          }
        return m.get(index);
        }

      inline bool has_entry(uint64_t index) const
        {
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == (uint64_t)-1)
            return false;
          if (v[ind].first == index)
            return true;
          }
        return m.has(index);
        }

    private:
      uint32_t mask;
      indexmap<T> m;
      std::vector<value_type> v;

      template <typename OtherT, typename sparse_buffer_type>
      friend class sparse_buffer_const_iterator;
    };


  ///////////////////////////////////////
  // concurrent_sparse_buffer
  ///////////////////////////////////////

  namespace concurrent_sparse_buffer_details
    {

    template <int A>
    struct get_power_two
      {
      static const int value = 2 * get_power_two<A - 1>::value;
      };

    template <>
    struct get_power_two<0>
      {
      static const int value = 1;
      };
    }


  template <typename T, typename concurrent_sparse_buffer_type>
  class concurrent_sparse_buffer_const_iterator
    {
    public:
      typedef concurrent_sparse_buffer_const_iterator<T, concurrent_sparse_buffer_type> self_type;
      typedef std::forward_iterator_tag iterator_category;
      typedef T value_type;
      typedef const T& reference;
      typedef std::ptrdiff_t difference_type;
      typedef typename std::pair<uint64_t, T> bucket;

      concurrent_sparse_buffer_const_iterator() : _si(nullptr) {}

      concurrent_sparse_buffer_const_iterator(const concurrent_sparse_buffer_type* si) : _si(si)
        {
        if (_si)
          {
          v_it = _si->v.begin();
          bucket_it = _si->table.begin();
          entry_it = bucket_it->begin();
          if (entry_it == bucket_it->end())
            _find_next_nonempty_bucket();
          }
        }

      reference operator*() const
        {
        return (v_it == _si->v.end()) ? *entry_it : v_it->second;
        }

      uint64_t entry() const
        {
        return (v_it == _si->v.end()) ? entry_it.key() : v_it->first;
        }

      self_type& operator++()
        {
        if (v_it == _si->v.end())
          {
          ++entry_it;
          if (entry_it == bucket_it->end())
            _find_next_nonempty_bucket();
          if (bucket_it == _si->table.end())
            {
            _si = nullptr;
            }
          }
        else
          {
          ++v_it;
          while (v_it != _si->v.end() && v_it->first == (uint64_t)-1)
            ++v_it;
          if (v_it == _si->v.end() && bucket_it == _si->table.end())
            _si = nullptr;
          }
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      bool operator == (const self_type& other) const
        {
        return ((_si == other._si) && (!_si || (v_it == other.v_it && bucket_it == other.bucket_it && entry_it == other.entry_it)));
        }

      bool operator != (const self_type& other) const
        {
        return !(*this == other);
        }

    private:

      void _find_next_nonempty_bucket()
        {
        ++bucket_it;
        while (bucket_it != _si->table.end() && bucket_it->empty())
          ++bucket_it;
        if (bucket_it != _si->table.end())
          {
          entry_it = bucket_it->begin();
          }
        }

    private:
      const concurrent_sparse_buffer_type* _si;
      typename std::vector<std::pair<uint64_t, T>>::const_iterator v_it;
      typename std::vector<flat_list_map<uint64_t, T>>::const_iterator bucket_it;
      typename flat_list_map<uint64_t, T>::const_iterator entry_it;
    };

  template <class T, int cluster_power = 1>
  class concurrent_sparse_buffer
    {
    public:
      typedef concurrent_sparse_buffer<T, cluster_power> self_type;
      typedef concurrent_sparse_buffer_const_iterator<T, concurrent_sparse_buffer<T, cluster_power>> const_iterator;
      typedef std::pair<uint64_t, T> value_type;
      typedef flat_list_map<uint64_t, T> bucket;

      concurrent_sparse_buffer(uint32_t nr_of_buckets = 4096)
        {
        if (nr_of_buckets == 0)
          nr_of_buckets = 4096;
        nr_of_buckets = get_nearest_power_of_two(nr_of_buckets);
        locks = new spinlock[nr_of_buckets];
        value_type dummy((uint64_t)-1, T());
        v.resize(nr_of_buckets << cluster_power, dummy);
        table.resize(nr_of_buckets);
        mask = nr_of_buckets - 1;
        }

      ~concurrent_sparse_buffer()
        {
        delete[] locks;
        }

      concurrent_sparse_buffer(const self_type& other) : v(other.v), mask(other.mask), table(other.table)
        {
        locks = new spinlock[mask + 1];
        }

      concurrent_sparse_buffer(self_type&& other) : concurrent_sparse_buffer()
        {
        swap(other);
        }

      self_type& operator=(const self_type& other)
        {
        self_type temp(other);
        swap(temp);
        return *this;
        }

      self_type& operator = (self_type&& other)
        {
            {
            self_type empty_temp;
            swap(empty_temp);
            }
            swap(other);
            return *this;
        }

      void swap(self_type& other)
        {
        std::swap(v, other.v);
        std::swap(table, other.table);
        std::swap(mask, other.mask);
        std::swap(locks, other.locks);
        }

      T& put(uint64_t index)
        {
        lock_bucket(index);
        T& res = get_entry(index);
        unlock_bucket(index);
        return res;
        }

      void erase(uint64_t index)
        {
        lock_bucket(index);
        /*
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < concurrent_sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == index)
            return v[ind].second;
          }
        return table[index & mask].get(index);
        */
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < concurrent_sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == index)
            {
            v[ind].first = (uint64_t)-1;
            for (int j = i + 1; j < concurrent_sparse_buffer_details::get_power_two<cluster_power>::value; ++j)
              {
              v[ind + j - 1] = v[ind + j];
              }
            unlock_bucket(index);
            return;
            }
          }
        table[index & mask].erase(index);
        unlock_bucket(index);
        }

      void put_safe(uint64_t index, const T& value)
        {
        lock_bucket(index);
        get_entry(index) = value;
        unlock_bucket(index);
        }

      const T& get(uint64_t index) const
        {
        lock_bucket(index);
        const T& res = get_entry(index);
        unlock_bucket(index);
        return res;
        }

      bool has(uint64_t index) const
        {
        lock_bucket(index);
        const bool res = has_entry(index);
        unlock_bucket(index);
        return res;
        }

      const_iterator begin() const
        {
        return const_iterator(this);
        }

      const_iterator end() const
        {
        return const_iterator(nullptr);
        }

      const_iterator cbegin() const
        {
        return const_iterator(this);
        }

      const_iterator cend() const
        {
        return const_iterator(nullptr);
        }

      uint64_t memory_used() const
        {
        uint64_t mem = table.size() * sizeof(bucket);
        for (const auto& b : table)
          {
          mem += b.memory_used();
          }
        mem += (mask + 1) * sizeof(spinlock);
        mem += v.capacity() * sizeof(value_type);
        return mem;
        }

      uint64_t size() const
        {
        uint64_t _sz = 0;
        for (const auto& pr : v)
          {
          if (pr.first != (uint64_t)(-1))
            ++_sz;
          }
        for (const auto& b : table)
          _sz += b.size();
        return _sz;
        }

      uint64_t bucket_count() const
        {
        return (uint64_t)(mask + 1);
        }

      double load_factor() const
        {
        return (double)size() / (double)bucket_count();
        }

    private:

      inline uint32_t get_nearest_power_of_two(uint32_t n)
        {
        --n;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        ++n;
        return n;
        }

      inline T& get_entry(uint64_t index)
        {
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < concurrent_sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == (uint64_t)-1)
            {
            v[ind].first = index;
            return v[ind].second;
            }
          if (v[ind].first == index)
            return v[ind].second;
          }
        return table[index & mask].put(index);
        }

      inline const T& get_entry(uint64_t index) const
        {
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < concurrent_sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == index)
            return v[ind].second;
          }
        return table[index & mask].get(index);
        }

      inline bool has_entry(uint64_t index) const
        {
        uint64_t ind = (index & mask) << cluster_power;
        for (int i = 0; i < concurrent_sparse_buffer_details::get_power_two<cluster_power>::value; ++i, ++ind)
          {
          if (v[ind].first == (uint64_t)-1)
            return false;
          if (v[ind].first == index)
            return true;
          }
        return table[index & mask].has(index);
        }

      inline void lock_bucket(uint64_t index) const
        {
        locks[index & mask].lock();
        }

      inline void unlock_bucket(uint64_t index) const
        {
        locks[index & mask].unlock();
        }

    private:
      uint32_t mask;
      std::vector<bucket> table;
      std::vector<value_type> v;
      mutable spinlock* locks;

      template <typename OtherT, typename concurrent_sparse_buffer_type>
      friend class concurrent_sparse_buffer_const_iterator;
    };


  //////////////////////////////////////////
  // circular_buffer
  //////////////////////////////////////////

  template <typename T, typename circular_buffer_type>
  class circular_buffer_iterator
    {
    public:

      typedef circular_buffer_iterator<T, circular_buffer_type> self_type;
      typedef std::random_access_iterator_tag     iterator_category;
      typedef typename circular_buffer_type::value_type      value_type;
      typedef typename circular_buffer_type::size_type       size_type;
      typedef typename circular_buffer_type::pointer         pointer;
      typedef typename circular_buffer_type::const_pointer   const_pointer;
      typedef typename circular_buffer_type::reference       reference;
      typedef typename circular_buffer_type::const_reference const_reference;
      typedef typename circular_buffer_type::difference_type difference_type;

      circular_buffer_iterator(circular_buffer_type *b, size_type p) : _buffer(b), _position(p) {}

      reference operator*() { return (*_buffer)[_position]; }
      pointer operator->() { return &(operator*()); }

      self_type& operator++()
        {
        _position += 1;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        _position -= 1;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      self_type operator + (difference_type n) const
        {
        self_type tmp(*this);
        tmp._position += n;
        return tmp;
        }

      self_type& operator+= (difference_type n)
        {
        _position += n;
        return *this;
        }

      self_type operator- (difference_type n) const
        {
        self_type tmp(*this);
        tmp._position -= n;
        return tmp;
        }

      self_type& operator-= (difference_type n)
        {
        _position -= n;
        return *this;
        }

      difference_type operator- (const self_type& c) const
        {
        return _position - c._position;
        }

      bool operator == (const self_type &other) const
        {
        return _position == other._position && _buffer == other._buffer;
        }

      bool operator != (const self_type& other) const
        {
        return _position != other._position && _buffer == other._buffer;
        }

      bool operator > (const self_type &other) const
        {
        return _position > other._position;
        }

      bool operator >= (const self_type &other) const
        {
        return _position >= other._position;
        }

      bool operator < (const self_type &other) const
        {
        return _position < other._position;
        }

      bool operator <= (const self_type &other) const
        {
        return _position <= other._position;
        }

    private:

      circular_buffer_type* _buffer;
      size_type  _position;
    };



  template <typename T, typename circular_buffer_type>
  class circular_buffer_const_iterator
    {
    public:

      typedef circular_buffer_const_iterator<T, circular_buffer_type> self_type;
      typedef std::random_access_iterator_tag     iterator_category;
      typedef typename circular_buffer_type::value_type      value_type;
      typedef typename circular_buffer_type::size_type       size_type;
      typedef typename circular_buffer_type::pointer         pointer;
      typedef typename circular_buffer_type::const_pointer   const_pointer;
      typedef typename circular_buffer_type::reference       reference;
      typedef typename circular_buffer_type::const_reference const_reference;
      typedef typename circular_buffer_type::difference_type difference_type;

      circular_buffer_const_iterator(const circular_buffer_type *b, size_type p) : _buffer(b), _position(p) {}

      const_reference operator*() { return (*_buffer)[_position]; }
      const_pointer operator->() { return &(operator*()); }

      self_type& operator++()
        {
        _position += 1;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        _position -= 1;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      self_type operator + (difference_type n) const
        {
        self_type tmp(*this);
        tmp._position += n;
        return tmp;
        }

      self_type& operator+= (difference_type n)
        {
        _position += n;
        return *this;
        }

      self_type operator- (difference_type n) const
        {
        self_type tmp(*this);
        tmp._position -= n;
        return tmp;
        }

      self_type& operator-= (difference_type n)
        {
        _position -= n;
        return *this;
        }

      difference_type operator- (const self_type& c) const
        {
        return _position - c._position;
        }

      bool operator == (const self_type &other) const
        {
        return _position == other._position && _buffer == other._buffer;
        }

      bool operator != (const self_type& other) const
        {
        return _position != other._position && _buffer == other._buffer;
        }

      bool operator > (const self_type &other) const
        {
        return _position > other._position;
        }

      bool operator >= (const self_type &other) const
        {
        return _position >= other._position;
        }

      bool operator < (const self_type &other) const
        {
        return _position < other._position;
        }

      bool operator <= (const self_type &other) const
        {
        return _position <= other._position;
        }

    private:

      const circular_buffer_type* _buffer;
      size_type  _position;
    };


  template <typename T, typename circular_buffer_type>
  class circular_buffer_const_reverse_iterator
    {
    public:

      typedef circular_buffer_const_reverse_iterator<T, circular_buffer_type> self_type;
      typedef std::random_access_iterator_tag     iterator_category;
      typedef typename circular_buffer_type::value_type      value_type;
      typedef typename circular_buffer_type::size_type       size_type;
      typedef typename circular_buffer_type::pointer         pointer;
      typedef typename circular_buffer_type::const_pointer   const_pointer;
      typedef typename circular_buffer_type::reference       reference;
      typedef typename circular_buffer_type::const_reference const_reference;
      typedef typename circular_buffer_type::difference_type difference_type;

      circular_buffer_const_reverse_iterator(const circular_buffer_type *b, size_type p) : _buffer(b), _position(p) {}

      const_reference operator*()
        {
        return (*_buffer)[_position - 1];
        }

      const_pointer operator->()
        {
        return &(operator*());
        }

      self_type& operator++()
        {
        _position -= 1;
        return *this;
        }

      self_type operator++(int)
        {
        self_type tmp(*this);
        ++(*this);
        return tmp;
        }

      self_type& operator--()
        {
        _position += 1;
        return *this;
        }

      self_type operator--(int)
        {
        self_type tmp(*this);
        --(*this);
        return tmp;
        }

      self_type operator + (difference_type n) const
        {
        self_type tmp(*this);
        tmp._position -= n;
        return tmp;
        }

      self_type& operator+= (difference_type n)
        {
        _position -= n;
        return *this;
        }

      self_type operator- (difference_type n) const
        {
        self_type tmp(*this);
        tmp._position += n;
        return tmp;
        }

      self_type& operator-= (difference_type n)
        {
        _position += n;
        return *this;
        }

      difference_type operator- (const self_type& c) const
        {
        return _position - c._position;
        }

      bool operator == (const self_type &other) const
        {
        return _position == other._position && _buffer == other._buffer;
        }

      bool operator != (const self_type& other) const
        {
        return _position != other._position && _buffer == other._buffer;
        }

      bool operator > (const self_type &other) const
        {
        return _position < other._position;
        }

      bool operator >= (const self_type &other) const
        {
        return _position <= other._position;
        }

      bool operator < (const self_type &other) const
        {
        return _position > other._position;
        }

      bool operator <= (const self_type &other) const
        {
        return _position >= other._position;
        }

    private:

      const circular_buffer_type* _buffer;
      size_type  _position;
    };

  template <typename T, typename alloc = std::allocator<T> >
  class circular_buffer
    {
    public:

      typedef circular_buffer<T, alloc>         self_type;

      typedef alloc                             allocator_type;

      typedef typename std::allocator_traits<alloc>::value_type          value_type;
      typedef typename std::allocator_traits<alloc>::pointer             pointer;
      typedef typename std::allocator_traits<alloc>::const_pointer       const_pointer;
      typedef typename std::allocator_traits<alloc>::value_type&         reference;
      typedef const typename std::allocator_traits<alloc>::value_type&   const_reference;

      typedef typename std::allocator_traits<alloc>::size_type         size_type;
      typedef typename std::allocator_traits<alloc>::difference_type   difference_type;

      typedef circular_buffer_iterator<T, self_type> iterator;
      typedef circular_buffer_const_iterator<T, self_type> const_iterator;
      typedef std::reverse_iterator<iterator>       reverse_iterator;
      typedef circular_buffer_const_reverse_iterator<T, self_type> const_reverse_iterator;

      explicit circular_buffer(size_type capacity = 1024, const alloc& allocator = alloc()) : _allocator(allocator), _array(_allocator.allocate(capacity)),
        _capacity(capacity), _head(1), _tail(0), _size(0)
        {
        }

      circular_buffer(const circular_buffer &other) : _array(_allocator.allocate(other._capacity)),
        _capacity(other._capacity), _head(other._head), _tail(other._tail), _size(other._size)
        {
        try
          {
          _Assign(other.begin(), other.end());
          }
        catch (...)
          {
          _Destroy();
          _allocator.deallocate(_array, _capacity);
          throw;
          }
        }

      template <class InputIterator>
      circular_buffer(InputIterator from, InputIterator to, const alloc& allocator = alloc()) : _allocator(allocator), _array(_allocator.allocate(std::distance(from, to))),
        _capacity(std::distance(from, to)), _head(1), _tail(0), _size(0)
        {
        circular_buffer tmp;
        tmp._Assign(from, to);
        swap(tmp);
        }

      ~circular_buffer()
        {
        _Destroy();
        _allocator.deallocate(_array, _capacity);
        }

      circular_buffer& operator = (const self_type& other)
        {
        circular_buffer tmp(other);
        swap(tmp);
        return *this;
        }

      void swap(circular_buffer& other)
        {
        std::swap(_array, other._array);
        std::swap(_capacity, other._capacity);
        std::swap(_head, other._head);
        std::swap(_tail, other._tail);
        std::swap(_size, other._size);
        }

      allocator_type get_allocator() const { return _allocator; }

      iterator         begin() { return iterator(this, 0); }
      iterator         end() { return iterator(this, size()); }

      const_iterator   begin() const { return const_iterator(this, 0); }
      const_iterator   end() const { return const_iterator(this, size()); }

      const_iterator   cbegin() const { return const_iterator(this, 0); }
      const_iterator   cend() const { return const_iterator(this, size()); }

      reverse_iterator rbegin() { return reverse_iterator(end()); }
      reverse_iterator rend() { return reverse_iterator(begin()); }

      const_reverse_iterator rbegin() const
        {
        return const_reverse_iterator(this, size());
        }

      const_reverse_iterator rend() const
        {
        return const_reverse_iterator(this, 0);
        }

      const_reverse_iterator crbegin() const
        {
        return const_reverse_iterator(this, size());
        }

      const_reverse_iterator crend() const
        {
        return const_reverse_iterator(this, 0);
        }

      size_type size() const { return _size; }

      size_type capacity() const { return _capacity; }

      bool empty() const { return !_size; }

      size_type max_size() const
        {
        return _allocator.max_size();
        }

      void reserve(size_type new_size)
        {
        if (capacity() < new_size)
          {
          circular_buffer tmp(new_size);
          tmp._Assign(begin(), end());
          swap(tmp);
          }
        }

      reference       front() { return _array[_head]; }
      reference       back() { return _array[_tail]; }
      const_reference front() const { return _array[_head]; }
      const_reference back() const { return _array[_tail]; }

      void push_back(const value_type &item)
        {
        size_type next = _Next_tail();
        if (_size == _capacity)
          {
          _array[next] = item;
          _Increment_head();
          }
        else
          {
          std::allocator_traits<alloc>::construct(_allocator, _array + next, item);
          }
        _Increment_tail();
        }

      void push_front(const value_type &item)
        {
        size_type previous = _Previous_head();
        if (_size == _capacity)
          {
          _array[previous] = item;
          _Decrement_tail();
          }
        else
          {
          std::allocator_traits<alloc>::construct(_allocator, _array + previous, item);
          }
        _Decrement_head();
        }

      void pop_back()
        {
        size_type destroy_pos = _tail;
        _Decrement_tail();
        std::allocator_traits<alloc>::destroy(_allocator, _array + destroy_pos);
        }

      void pop_front()
        {
        size_type destroy_pos = _head;
        _Increment_head();
        std::allocator_traits<alloc>::destroy(_allocator, _array + destroy_pos);
        }

      void clear()
        {
        for (size_type i = 0; i < _size; ++i)
          {
          std::allocator_traits<alloc>::destroy(_allocator, _array + _Index(i));
          }
        _head = 1;
        _tail = 0;
        _size = 0;
        }

      reference       operator [] (size_type i) { return _array[_Index(i)]; }
      const_reference operator [] (size_type i) const { return _array[_Index(i)]; }

      reference       at(size_type i) { return _At(i); }
      const_reference at(size_type i) const { return _At(i); }

    private:

      reference _At(size_type index) const
        {
        if (index >= _size)
          {
          throw std::out_of_range("circular_buffer indexed access is out of range.");
          }
        return _array[_Index(index)];
        }

      size_type _Index(size_type index) const
        {
        return (index + _head) % _capacity;
        }

      void _Increment_tail()
        {
        ++_size;
        _tail = _Next_tail();
        }

      void _Decrement_tail()
        {
        --_size;
        _tail = _Previous_tail();
        }

      size_type _Next_tail()
        {
        return (_tail + 1 == _capacity) ? 0 : _tail + 1;
        }

      size_type _Previous_tail()
        {
        return (_tail == 0) ? _capacity - 1 : _tail - 1;
        }

      size_type _Previous_head()
        {
        return (_head == 0) ? _capacity - 1 : _head - 1;
        }

      void _Increment_head()
        {
        ++_head;
        --_size;
        if (_head == _capacity)
          _head = 0;
        }

      void _Decrement_head()
        {
        ++_size;
        _head = _Previous_head();
        }

      template <typename f_iter>
      void _Assign(f_iter from, f_iter to)
        {
        if (_size) clear();
        while (from != to)
          {
          push_back(*from);
          ++from;
          }
        }

      void _Destroy()
        {
        for (size_type i = 0; i < _size; ++i)
          {
          std::allocator_traits<alloc>::destroy(_allocator, _array + _Index(i));
          }
        }

      allocator_type  _allocator;
      value_type*     _array;
      size_type       _capacity;
      size_type       _head;
      size_type       _tail;
      size_type       _size;
    };

  template <typename T, typename alloc>
  bool operator == (const circular_buffer<T, alloc>& a, const circular_buffer<T, alloc> &b)
    {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
    }

  template <typename T, typename alloc>
  bool operator != (const circular_buffer<T, alloc>& a, const circular_buffer<T, alloc>& b)
    {
    return a.size() != b.size() || !std::equal(a.begin(), a.end(), b.begin());
    }

  template <typename T, typename alloc>
  bool operator < (const circular_buffer<T, alloc>& a, const circular_buffer<T, alloc>& b)
    {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }

  template <typename T, typename alloc>
  bool operator <= (const circular_buffer<T, alloc>& a, const circular_buffer<T, alloc>& b)
    {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end(), [](const T& left, const T& right) {return left <= right; });
    }

  template <typename T, typename alloc>
  bool operator > (const circular_buffer<T, alloc>& a, const circular_buffer<T, alloc>& b)
    {
    return !(b <= a);
    }

  template <typename T, typename alloc>
  bool operator >= (const circular_buffer<T, alloc>& a, const circular_buffer<T, alloc>& b)
    {
    return !(b < a);
    }

  /////////////////////////////////////////////////
  // circular_queue
  /////////////////////////////////////////////////

  template <typename T>
  class circular_queue
    {
    public:
      typedef circular_queue<T> self_type;
      typedef circular_buffer<T> container_type;
      typedef typename container_type::value_type value_type;
      typedef typename container_type::size_type size_type;
      typedef typename container_type::reference reference;
      typedef typename container_type::const_reference const_reference;

      circular_queue(size_type capacity = 1024) : c(capacity)
        {
        }

      circular_queue(const self_type& _Right) : c(_Right.c)
        {
        }

      explicit circular_queue(const container_type& _Cont) : c(_Cont)
        {
        }

      self_type& operator = (const self_type& _Right)
        {
        c = _Right.c;
        return (*this);
        }

      void swap(self_type& _Right)
        {
        c.swap(_Right.c);
        }

      circular_queue(self_type&& _Right) : c(std::move(_Right.c))
        {
        }

      explicit circular_queue(container_type&& _Cont) : c(std::move(_Cont))
        {
        }

      self_type& operator= (self_type&& _Right)
        {
        c = std::move(_Right.c);
        return (*this);
        }

      void push(value_type&& _Val)
        {
        c.push_back(std::move(_Val));
        }

      bool empty() const
        {
        return (c.empty());
        }

      size_type size() const
        {
        return (c.size());
        }

      size_type capacity() const { return c.capacity(); }

      void reserve(size_type new_size)
        {
        c.reserve(new_size);
        }

      reference front()
        {
        return (c.front());
        }

      const_reference front() const
        {
        return (c.front());
        }

      reference back()
        {
        return (c.back());
        }

      const_reference back() const
        {
        return (c.back());
        }

      void push(const value_type& _Val)
        {
        c.push_back(_Val);
        }

      void pop()
        {
        c.pop_front();
        }

      void clear()
        {
        c.clear();
        }

    private:
      container_type c;

    };



  } // namespace jtk
