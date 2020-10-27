#include "container_tests.h"
#include "test_assert.h"

#include "../jtk/containers.h"
#include <array>

#include <iostream>

namespace jtk
  {

  void TEST_HashedHeapBasicTest()
    {
    std::vector<std::pair<size_t, double>> vec;
    vec.push_back({ 0, 1. });
    vec.push_back({ 1, 2. });
    vec.push_back({ 4, -3. });
    vec.push_back({ 10, 1.5 });
    vec.push_back({ 3, 4. });

    hashed_heap<size_t, double> heap;
    for (auto v : vec)
      heap.push(v);

    std::sort(vec.begin(), vec.end(), [](auto p0, auto p1) { return p0.second < p1.second; });
    size_t cnt = 0;
    while (!heap.empty())
      {
      auto t = heap.top();
      heap.pop();
      TEST_ASSERT(t == vec[cnt++]);
      }
    }

  void TEST_HashedHeapRepeatedLessTest()
    {
    std::vector<std::pair<size_t, double>> vec;
    vec.push_back({ 0, 1. });
    vec.push_back({ 1, 2. });
    vec.push_back({ 4, -3. });
    vec.push_back({ 10, 1.5 });
    vec.push_back({ 3, 4. });
    vec.push_back({ 0, .4 });

    hashed_heap<size_t, double> heap;
    for (auto v : vec)
      heap.push(v);

    vec.erase(vec.begin());
    std::sort(vec.begin(), vec.end(), [](auto p0, auto p1) { return p0.second < p1.second; });
    size_t cnt = 0;
    while (!heap.empty())
      {
      auto t = heap.top();
      heap.pop();
      TEST_ASSERT(t == vec[cnt++]);
      }
    }

  void TEST_HashedHeapRepeatedGreaterTest()
    {
    std::vector<std::pair<size_t, double>> vec;
    vec.push_back({ 0, .4 });
    vec.push_back({ 1, 2. });
    vec.push_back({ 4, -3. });
    vec.push_back({ 10, 1.5 });
    vec.push_back({ 3, 4. });
    vec.push_back({ 0, 1. });

    hashed_heap<size_t, double> heap;
    for (auto v : vec)
      heap.push(v);

    vec.erase(vec.end() - 1);
    std::sort(vec.begin(), vec.end(), [](auto p0, auto p1) { return p0.second < p1.second; });
    size_t cnt = 0;
    while (!heap.empty())
      {
      auto t = heap.top();
      heap.pop();
      TEST_ASSERT(t == vec[cnt++]);
      }
    }

  void TEST_HashedHeapNonDefaultKeyTest()
    {
    typedef std::array<size_t, 2> edge;
    struct edge_equal
      {
      bool operator() (const edge& e0, const edge& e1) const
        {
        return e0 == e1 || (e0[0] == e1[1] && e0[1] == e1[0]);
        }
      };
    struct edge_hash
      {
      size_t operator()(const edge& e) const
        {
        return (e[0] + e[1]) * (e[0] + e[1] + 1) / 2;
        }
      };

    std::vector<std::pair<edge, double>> vec;
    vec.push_back({ { { 0, 1 } }, .4 });
    vec.push_back({ { { 1, 4 } }, 2. });
    vec.push_back({ { { 4, 3 } }, -3. });
    vec.push_back({ { { 10, 15 } }, 1.5 });
    vec.push_back({ { { 3, 7 } }, 4. });
    vec.push_back({ { { 1, 0 } }, 1. });

    hashed_heap<edge, double, edge_hash, edge_equal> heap;
    for (auto v : vec)
      heap.push(v);

    vec.erase(vec.end() - 1);
    std::sort(vec.begin(), vec.end(), [](auto p0, auto p1) { return p0.second < p1.second; });
    size_t cnt = 0;
    while (!heap.empty())
      {
      auto t = heap.top();
      heap.pop();
      TEST_ASSERT(t == vec[cnt++]);
      }
    }

  void TEST_HashedHeapNonDefaultKeyAndDataTest()
    {
    typedef std::array<size_t, 2> edge;
    struct edge_equal
      {
      bool operator() (const edge& e0, const edge& e1) const
        {
        return e0 == e1 || (e0[0] == e1[1] && e0[1] == e1[0]);
        }
      };
    struct edge_hash
      {
      size_t operator()(const edge& e) const
        {
        return (e[0] + e[1]) * (e[0] + e[1] + 1) / 2;
        }
      };

    typedef std::array<double, 2> data;
    struct data_less
      {
      bool operator() (const data& d0, const data& d1) const
        {
        return d0[1] < d1[1];
        }
      };

    std::vector<std::pair<edge, data>> vec;
    vec.push_back({ { { 0, 1 } }, { { .4, .2 } } });
    vec.push_back({ { { 1, 4 } }, { { 2., 1. } } });
    vec.push_back({ { { 4, 3 } }, { { -3., 5. } } });
    vec.push_back({ { { 10, 15 } }, { { 1.5, -4.5 } } });
    vec.push_back({ { { 3, 7 } }, { { 4., -9. } } });
    vec.push_back({ { { 1, 0 } }, { { 1., 3. } } });

    hashed_heap < edge, data, edge_hash, edge_equal, data_less > heap;
    for (auto v : vec)
      heap.push(v);

    vec.erase(vec.end() - 1);
    std::sort(vec.begin(), vec.end(), [](auto p0, auto p1) { return p0.second[1] < p1.second[1]; });
    size_t cnt = 0;
    while (!heap.empty())
      {
      auto t = heap.top();
      heap.pop();
      TEST_ASSERT(t == vec[cnt++]);
      }
    }


  void TEST_aligned_vector_test_1()
    {
    aligned_vector<double> v(5);
    TEST_EQ(5, v.size());
    v.push_back(3.0);
    TEST_EQ(3.0, v[5]);
    aligned_vector<double> v2;
    TEST_ASSERT(v2.empty());
    }

  void TEST_aligned_vector_test_2()
    {
    aligned_vector<double> v;
    v.reserve(100);
    TEST_EQ(100, v.capacity());
    }

  void TEST_aligned_vector_test_3()
    {
    aligned_vector<double> v;
    v.push_back(1.0);
    v.push_back(2.0);
    v.push_back(3.0);
    v.push_back(4.0);
    v.push_back(5.0);
    auto it = v.begin();
    auto it_end = v.end();
    double val = 1.0;
    for (; it != it_end; ++it)
      {
      TEST_EQ(val, *it);
      val += 1.0;
      }

    aligned_vector<double> v2(v);
    it = v.begin();
    val = 1.0;
    for (; it != it_end; ++it)
      {
      TEST_EQ(val, *it);
      val += 1.0;
      }

    for (auto& value : v2)
      value *= 2.0;
    it = v2.begin();
    it_end = v2.end();
    val = 1.0;
    for (; it != it_end; ++it)
      {
      TEST_EQ(val*2.0, *it);
      val += 1.0;
      }
    }

  void TEST_aligned_vector_test_4()
    {
    aligned_vector<double> v(5, 3.14);
    TEST_EQ(5, v.size());
    TEST_EQ(3.14, v[0]);
    TEST_EQ(3.14, v[1]);
    TEST_EQ(3.14, v[2]);
    TEST_EQ(3.14, v[3]);
    TEST_EQ(3.14, v[4]);
    }

  void TEST_typed_memory_pool_init()
    {
    typed_memory_pool<uint64_t, 16> pool;
    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    uint64_t* entry = pool.allocate();
    *entry = 54;

    TEST_EQ((1 << 16) - 1, pool.available());
    TEST_EQ(1, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    pool.deallocate(entry);

    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    entry = pool.allocate();
    TEST_EQ(54, *entry);
    }

  void TEST_concurrent_typed_memory_pool_init()
    {
    concurrent_typed_memory_pool<uint64_t, 16> pool;
    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    uint64_t* entry = pool.allocate();
    *entry = 54;

    TEST_EQ((1 << 16) - 1, pool.available());
    TEST_EQ(1, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    pool.deallocate(entry);

    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    entry = pool.allocate();
    TEST_EQ(54, *entry);
    }


  void TEST_memory_pool_init()
    {
    memory_pool<sizeof(uint64_t), 16> pool;
    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    uint8_t* entry = pool.allocate();
    const uint64_t value = 12345678987654;
    memcpy(entry, &value, sizeof(uint64_t));

    TEST_EQ((1 << 16) - 1, pool.available());
    TEST_EQ(1, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    pool.deallocate(entry);

    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    entry = pool.allocate();
    uint64_t actual;
    memcpy(&actual, entry, sizeof(uint64_t));
    TEST_EQ(value, actual);
    }

  void TEST_concurrent_memory_pool_init()
    {
    concurrent_memory_pool<sizeof(uint64_t), 16> pool;
    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    uint8_t* entry = pool.allocate();
    const uint64_t value = 12345678987654;
    memcpy(entry, &value, sizeof(uint64_t));

    TEST_EQ((1 << 16) - 1, pool.available());
    TEST_EQ(1, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    pool.deallocate(entry);

    TEST_EQ(1 << 16, pool.available());
    TEST_EQ(0, pool.allocated());
    TEST_EQ(1 << 16, pool.size());

    entry = pool.allocate();
    uint64_t actual;
    memcpy(&actual, entry, sizeof(uint64_t));
    TEST_EQ(value, actual);
    }

  void flat_map_1()
    {
    flat_map<uint64_t, uint64_t> m;

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    for (int i = 0; i < 5000; ++i)
      TEST_EQ(i * 3, m.get(i));
    }

  void flat_map_2()
    {
    flat_map<uint64_t, uint64_t> m;

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    auto cit = m.cbegin();
    auto cend = m.cend();

    for (; cit != cend; ++cit)
      {
      TEST_EQ(cit.key() * 3, *cit);
      }

    TEST_ASSERT(m.size() == 5000);
    }

  void flat_map_3()
    {
    flat_map<uint64_t, uint64_t> m;

    for (int i = 0; i < 5000; ++i)
      {
      if (i % 2)
        {
        m.put(i) = i * 3;
        }
      }

    auto cit = m.begin();
    auto cend = m.end();

    for (; cit != cend; ++cit)
      {
      *cit *= 2;
      }

    cit = m.begin();
    cend = m.end();

    for (; cit != cend; ++cit)
      {
      TEST_EQ(cit.key() * 6, *cit);
      }

    TEST_ASSERT(m.size() == 2500);
    }

  void hash_1()
    {
    indexmap<uint64_t> m;

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    for (int i = 0; i < 5000; ++i)
      TEST_EQ(i * 3, m.get(i));
    }

  void hash_2()
    {
    indexmap<uint64_t> m;

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    auto cit = m.cbegin();
    auto cend = m.cend();

    for (; cit != cend; ++cit)
      {
      TEST_EQ(cit.entry() * 3, *cit);
      }
    }

  void sparse_vector_1()
    {
    sparse_vector<uint64_t> m(1024);

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    for (int i = 0; i < 5000; ++i)
      TEST_EQ(i * 3, m.get(i));
    }

  void sparse_vector_2()
    {
    sparse_vector<uint64_t> m(1024);

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    auto cit = m.cbegin();
    auto cend = m.cend();

    for (; cit != cend; ++cit)
      {
      TEST_EQ(cit.entry() * 3, *cit);
      }
    }

  void concurrent_sparse_vector_1()
    {
    concurrent_sparse_vector<uint64_t> m(1024);

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    for (int i = 0; i < 5000; ++i)
      TEST_EQ(i * 3, m.get(i));
    }

  void concurrent_sparse_vector_2()
    {
    concurrent_sparse_vector<uint64_t> m(1024);

    for (int i = 0; i < 5000; ++i)
      m.put(i) = i * 3;

    auto cit = m.cbegin();
    auto cend = m.cend();

    for (; cit != cend; ++cit)
      {
      TEST_EQ(cit.entry() * 3, *cit);
      }
    }


  void test_circular_buffer_construction()
    {
    circular_buffer<int> cbuf(20);
    TEST_EQ(0, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(cbuf.empty());
    }

  void test_circular_buffer_pushback()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    }

  void test_circular_buffer_popback()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_back();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    }

  void test_circular_buffer_at()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(cbuf.at(i), i * 2);
    try
      {
      std::cout << cbuf.at(50) << "This should not print\n";
      }
    catch (const std::out_of_range& e)
      {
      TEST_ASSERT(std::string(e.what()) == std::string("circular_buffer indexed access is out of range."));
      }
    }

  void test_circular_buffer_pop_front()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    }

  void test_circular_buffer_adding_over_capacity_after_pop_front()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);
    }


  void test_circular_buffer_push_front()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 20; ++i)
      cbuf.push_front(i * 2);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 20; ++i)
      TEST_EQ((19 - i) * 2, cbuf[i]);

    for (int i = 0; i < 5; ++i)
      cbuf.push_front((i + 20) * 2);

    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 20) * 2, cbuf[4 - i]);

    for (int i = 5; i < 20; ++i)
      TEST_EQ((19 - i + 5) * 2, cbuf[i]);
    }

  void test_circular_buffer_pop_back_over_capacity()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    for (int i = 0; i < 5; ++i)
      cbuf.pop_back();

    TEST_EQ(15, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 15; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);
    }

  void test_circular_buffer_iterator()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);


    circular_buffer<int>::iterator it = cbuf.begin();
    int i = 0;
    while (it != cbuf.end())
      {
      if (i < 4)
        TEST_EQ((i + 6) * 2, *it);
      else
        TEST_EQ((i - 4) * 2 + 1, *it);
      ++it;
      ++i;
      }

    circular_buffer<int>::const_iterator cit = cbuf.cbegin();
    i = 0;
    while (cit != cbuf.cend())
      {
      if (i < 4)
        TEST_EQ((i + 6) * 2, *cit);
      else
        TEST_EQ((i - 4) * 2 + 1, *cit);
      ++cit;
      ++i;
      }
    }


  void test_circular_buffer_const_iterator()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    const circular_buffer<int>& const_cbuf(cbuf);
    circular_buffer<int>::const_iterator it = const_cbuf.begin();
    int i = 0;
    while (it != const_cbuf.end())
      {
      if (i < 4)
        TEST_EQ((i + 6) * 2, *it);
      else
        TEST_EQ((i - 4) * 2 + 1, *it);
      ++it;
      ++i;
      }

    }


  void test_circular_buffer_reverse_iterator()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int>::reverse_iterator it = cbuf.rbegin();
    int i = 19;
    while (it != cbuf.rend())
      {
      if (i < 4)
        TEST_EQ((i + 6) * 2, *it);
      else
        TEST_EQ((i - 4) * 2 + 1, *it);
      ++it;
      --i;
      }

    circular_buffer<int>::const_reverse_iterator cit = cbuf.crbegin();
    i = 19;
    while (cit != cbuf.crend())
      {
      if (i < 4)
        TEST_EQ((i + 6) * 2, *cit);
      else
        TEST_EQ((i - 4) * 2 + 1, *cit);
      ++cit;
      --i;
      }
    }


  void test_circular_buffer_copy_assignment()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int> cbuf2 = cbuf;
    TEST_ASSERT(cbuf == cbuf2);
    TEST_ASSERT(!(cbuf != cbuf2));
    cbuf2.push_back(0);
    TEST_ASSERT(!(cbuf == cbuf2));
    TEST_ASSERT(cbuf != cbuf2);
    }


  void test_circular_buffer_copy_constructor()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int> cbuf2(cbuf);
    TEST_ASSERT(cbuf == cbuf2);
    TEST_ASSERT(!(cbuf != cbuf2));
    cbuf2.push_back(0);
    TEST_ASSERT(!(cbuf == cbuf2));
    TEST_ASSERT(cbuf != cbuf2);
    }


  void test_circular_buffer_assignment()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int> cbuf2;
    cbuf2 = cbuf;
    TEST_ASSERT(cbuf == cbuf2);
    TEST_ASSERT(!(cbuf != cbuf2));
    cbuf2.push_back(0);
    TEST_ASSERT(!(cbuf == cbuf2));
    TEST_ASSERT(cbuf != cbuf2);
    }


  void test_circular_buffer_reserve()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int> cbuf2;
    cbuf2 = cbuf;
    cbuf2.reserve(1024);
    TEST_EQ(1024, cbuf2.capacity());
    TEST_ASSERT(cbuf == cbuf2);
    TEST_ASSERT(!(cbuf != cbuf2));
    cbuf2.push_back(0);
    TEST_ASSERT(!(cbuf == cbuf2));
    TEST_ASSERT(cbuf != cbuf2);
    }


  void test_circular_buffer_ranged_construction()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int> cbuf2(cbuf.begin(), cbuf.end());
    TEST_ASSERT(cbuf == cbuf2);
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf2[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf2[i]);
    }

  void test_circular_buffer_clear()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);
    cbuf.clear();
    TEST_EQ(0, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(cbuf.empty());
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2 + 1, cbuf[i]);
    }


  void test_circular_buffer_swap()
    {
    circular_buffer<int> cbuf(20);
    for (int i = 0; i < 10; ++i)
      cbuf.push_back(i * 2);
    TEST_EQ(10, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 10; ++i)
      TEST_EQ(i * 2, cbuf[i]);
    for (int i = 0; i < 5; ++i)
      cbuf.pop_front();
    TEST_EQ(5, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 5; ++i)
      TEST_EQ((i + 5) * 2, cbuf[i]);
    for (int i = 0; i < 16; ++i)
      cbuf.push_back(i * 2 + 1);
    TEST_EQ(20, cbuf.size());
    TEST_EQ(20, cbuf.capacity());
    TEST_ASSERT(!cbuf.empty());
    for (int i = 0; i < 4; ++i)
      TEST_EQ((i + 6) * 2, cbuf[i]);
    for (int i = 4; i < 20; ++i)
      TEST_EQ((i - 4) * 2 + 1, cbuf[i]);

    circular_buffer<int> cbuf2(0);
    std::swap(cbuf2, cbuf);

    TEST_EQ(0, cbuf.size());
    TEST_EQ(0, cbuf.capacity());
    TEST_ASSERT(cbuf.empty());

    cbuf.reserve(10);

    TEST_EQ(0, cbuf.size());
    TEST_EQ(10, cbuf.capacity());
    TEST_ASSERT(cbuf.empty());
    }

  void test_circular_queue_construction()
    {
    circular_queue<int> cq(20);
    TEST_EQ(0, cq.size());
    TEST_EQ(20, cq.capacity());
    TEST_ASSERT(cq.empty());
    }

  void test_circular_queue_push()
    {
    circular_queue<int> cq(20);
    for (int i = 0; i < 10; ++i)
      cq.push(i * 2);
    for (int i = 0; i < 10; ++i)
      {
      TEST_EQ(i * 2, cq.front());
      cq.pop();
      }
    }

  void test_circular_queue_push_over_limit()
    {
    circular_queue<int> cq(20);
    for (int i = 0; i < 30; ++i)
      cq.push(i * 2);
    for (int i = 0; i < 20; ++i)
      {
      TEST_EQ((i + 10) * 2, cq.front());
      cq.pop();
      }
    }

  } // namespace jtk

void run_all_container_tests()
  {
  using namespace jtk;
  TEST_HashedHeapBasicTest();
  TEST_HashedHeapRepeatedLessTest();
  TEST_HashedHeapRepeatedGreaterTest();
  TEST_HashedHeapNonDefaultKeyTest();
  TEST_HashedHeapNonDefaultKeyAndDataTest();

  TEST_aligned_vector_test_1();
  TEST_aligned_vector_test_2();
  TEST_aligned_vector_test_3();
  TEST_aligned_vector_test_4();

  TEST_typed_memory_pool_init();
  TEST_concurrent_typed_memory_pool_init();
  TEST_memory_pool_init();
  TEST_concurrent_memory_pool_init();

  flat_map_1();
  flat_map_2();
  flat_map_3();

  hash_1();
  hash_2();
  sparse_vector_1();
  sparse_vector_2();
  concurrent_sparse_vector_1();
  concurrent_sparse_vector_2();

  test_circular_buffer_construction();
  test_circular_buffer_pushback();
  test_circular_buffer_popback();
  test_circular_buffer_at();
  test_circular_buffer_pop_front();
  test_circular_buffer_adding_over_capacity_after_pop_front();
  test_circular_buffer_push_front();
  test_circular_buffer_pop_back_over_capacity();
  test_circular_buffer_iterator();
  test_circular_buffer_const_iterator();
  test_circular_buffer_reverse_iterator();
  test_circular_buffer_copy_assignment();
  test_circular_buffer_copy_constructor();
  test_circular_buffer_assignment();
  test_circular_buffer_reserve();
  test_circular_buffer_ranged_construction();
  test_circular_buffer_clear();
  test_circular_buffer_swap();

  test_circular_queue_construction();
  test_circular_queue_push();
  test_circular_queue_push_over_limit();
  }