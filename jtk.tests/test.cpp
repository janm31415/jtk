#include "test_assert.h"
#include "concurrency_tests.h"
#include "container_tests.h"
#include "file_utils_tests.h"
#include "geometry_tests.h"
#include "image_tests.h"
#include "log_tests.h"
#include "mat_tests.h"
#include "qbvh_tests.h"
#include "vec_tests.h"

#include <ctime>

int main(int /*argc*/, const char* /*argv*/[])
  {
  InitTestEngine();

  auto tic = std::clock();
  run_all_geometry_tests();
  run_all_concurrency_tests();
  run_all_container_tests();
  run_all_file_utils_tests();
  run_all_image_tests();
  run_all_mat_tests();
  run_all_qbvh_tests();
  run_all_vec_tests();  
  run_all_log_tests();
  auto toc = std::clock();

  if (!testing_fails) 
    {
    TEST_OUTPUT_LINE("Succes: %d tests passed.", testing_success);
    }
  else 
    {
    TEST_OUTPUT_LINE("FAILURE: %d out of %d tests failed (%d failures).", testing_fails, testing_success+testing_fails, testing_fails);
    }
  TEST_OUTPUT_LINE("Test time: %f seconds.", (double)(toc - tic)/(double)CLOCKS_PER_SEC);
  return CloseTestEngine(true);
  }
