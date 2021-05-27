#include "log_tests.h"
#include "test_assert.h"

#define JTK_LOG_IMPLEMENTATION
#include "../jtk/log.h"

#include <iostream>
#include <sstream>

namespace jtk
  {

  void test_time()
    {
    std::string str = now_time();
    TEST_ASSERT(!str.empty());
    }

  void test_log()
    {    
    std::stringstream ss;
    init_log_stream(&ss, false);
    log().get(log_level::info) << "Test";
    TEST_EQ(std::string(" INFO: Test\n"), ss.str());
    release_log_stream();
    }

  }

void run_all_log_tests()
  {
  using namespace jtk;
  test_time();
  test_log();
  }