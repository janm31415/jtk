#include "concurrency_tests.h"

#include "../jtk/concurrency.h"
#include <thread>
#include <stdint.h>

#include "test_assert.h"

#include <iostream>

namespace jtk
  {

  spinlock_rw g_lock;
  uint64_t g_counter = 0;
  const uint32_t count = 5000;

  void add_job()
    {
    for (int i = 0; i < count; ++i)
      {
      g_lock.lock();
      ++g_counter;
      std::cout << "Thread " << get_thread_id() << std::endl;
      g_lock.unlock();
      }
    }

  void read_job()
    {
    for (int i = 0; i < count; ++i)
      {
      g_lock.lock_read();
      std::cout << g_counter << std::endl;
      g_lock.unlock();
      }
    }

  void spin_lock_rw_test()
    {
    std::thread th1(add_job);
    std::thread th2(add_job);

    std::thread th3(read_job);
    std::thread th4(read_job);
    std::thread th5(read_job);

    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();

    TEST_EQ(10000, g_counter);
    }

  }

void run_all_concurrency_tests()
  {
  using namespace jtk;
  spin_lock_rw_test();
  }