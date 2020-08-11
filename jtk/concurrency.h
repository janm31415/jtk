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
#elif defined(_ENABLE_THREADS)
#include <thread>
#include <cassert>
#include <exception>
#include <vector>
#else
#include <vector>
#include <cassert>
#include <exception>
#endif

namespace jtk
  {

  inline unsigned int hardware_concurrency()
    {
#if defined(_ENABLE_TBB)
    return std::thread::hardware_concurrency();
#elif defined(_ENABLE_PPL)
    return std::thread::hardware_concurrency();
#elif defined(_ENABLE_THREADS)
    return std::thread::hardware_concurrency();
#else
    return 1;
#endif
    }


  template <class _Type, class TFunctor>
  void parallel_for(_Type first, _Type last, TFunctor fun)
    {
#if defined(_ENABLE_TBB)
    tbb::parallel_for(first, last, fun);
#elif defined(_ENABLE_PPL)
    Concurrency::parallel_for(first, last, fun);
#elif defined(_ENABLE_THREADS)

    const _Type n_threads = (_Type)std::thread::hardware_concurrency();

    const _Type n = last - first;

    const _Type n_max_tasks_per_thread = (n / n_threads) + (n % n_threads == 0 ? 0 : 1);
    const _Type n_lacking_tasks = n_max_tasks_per_thread * n_threads - n;

    auto inner_loop = [&](const _Type thread_index)
      {
      const _Type n_lacking_tasks_so_far = n_threads > (thread_index + n_lacking_tasks) ? 0 : (thread_index + n_lacking_tasks) - n_threads;
      const _Type inclusive_start_index = thread_index * n_max_tasks_per_thread - n_lacking_tasks_so_far;
      const _Type exclusive_end_index = inclusive_start_index + n_max_tasks_per_thread - (thread_index + n_lacking_tasks >= n_threads ? 1 : 0);

      for (_Type k = inclusive_start_index; k < exclusive_end_index; ++k)
        {
        fun(k + first);
        }
      };
    std::vector<std::thread> threads;
    for (_Type j = 0; j < n_threads; ++j) { threads.push_back(std::thread(inner_loop, j)); }
    for (auto& t : threads) { t.join(); }
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
#elif defined(_ENABLE_THREADS)
    std::sort(begin, end, comp);
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
#elif defined(_ENABLE_THREADS)
    std::sort(begin, end);
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
#elif defined(_ENABLE_THREADS)
  // no concurrent_vector implementation
#else
  template<typename T>
  using concurrent_vector = std::vector<T>;
#endif

  } // namespace jtk