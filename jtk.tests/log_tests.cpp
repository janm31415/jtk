#include "log_tests.h"
#include "test_assert.h"

#define JTK_LOG_IMPLEMENTATION
#include "../jtk/log.h"

#include <iostream>

namespace jtk
  {

  void test_time()
    {
    std::string str = now_time();
    TEST_ASSERT(!str.empty());
    }

  void test_log()
    {
    init_log_stream(&std::cout);
    log().get(log_level::info) << "Test";
    }

  }

void run_all_log_tests()
  {
  using namespace jtk;
  test_time();
  test_log();
  }