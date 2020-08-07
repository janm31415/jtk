#pragma once

#if defined(_ENABLE_TBB)
#include <tbb/enumerable_thread_specific.h>
#undef min
#undef max
#include <thread>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/concurrent_vector.h>
#include <tbb/spin_rw_mutex.h>
#include <cassert>
#include <exception>
#elif defined(_ENABLE_PPL)
#include <ppl.h>
#include <thread>
#include <cassert>
#include <exception>
#include <ppl.h>
#include <concurrent_vector.h>
#else
#include <vector>
#include <cassert>
#include <exception>
#endif;

namespace jtk
  {

  inline unsigned int hardware_concurrency()
    {
#if defined(_ENABLE_TBB)
    return std::thread::hardware_concurrency();
    #elif defined(_ENABLE_PPL)
      return std::thread::hardware_concurrency();
#else
    return 1;
#endif;
    }


  template <class _Type, class TFunctor>
  void parallel_for(_Type first, _Type last, TFunctor fun)
    {
#if defined(_ENABLE_TBB)
    tbb::parallel_for(first, last, fun);
#elif defined(_ENABLE_PPL)
    Concurrency::parallel_for(first, last, fun);
#else
    for (_Type i = first; i != last; ++i)
      fun(i);
#endif
    }

  template<typename RandomAccessIterator, typename Compare>
  inline void parallel_sort(RandomAccessIterator begin, RandomAccessIterator end, const Compare& comp)
    {
#if defined(_ENABLE_TBB)
    tbb::parallel_sort(begin, end, comp);
#elif defined(_ENABLE_PPL)
    Concurrency::parallel_sort(begin, end, comp);
#else
    std::sort(begin, end, comp);
#endif
    }

  template<typename RandomAccessIterator>
  inline void parallel_sort(RandomAccessIterator begin, RandomAccessIterator end)
    {
#if defined(_ENABLE_TBB)
    tbb::parallel_sort(begin, end);
#elif defined(_ENABLE_PPL)
    Concurrency::parallel_sort(begin, end);
#else
    std::sort(begin, end);
#endif
    }


#if defined(_ENABLE_TBB)
  template<typename T>
  using concurrent_vector = tbb::concurrent_vector<T>;
#elif defined(_ENABLE_PPL)
  template<typename T>
  using concurrent_vector = Concurrency::concurrent_vector<T>;
#else
  template<typename T>
  using concurrent_vector = std::vector<T>;
#endif

  } // namespace jtk