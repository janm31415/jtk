#include "file_utils_tests.h"
#include "test_assert.h"

#include "../jtk/file_utils.h"

namespace jtk
  {

  void csv_write_and_read()
    {
    std::vector<std::vector<std::string>> data;
    std::vector<std::string> line1, line2;
    line1.emplace_back("The");
    line1.emplace_back("first");
    line1.emplace_back("line");
    line2.emplace_back("And");
    line2.emplace_back("the");
    line2.emplace_back("second");
    line2.emplace_back("line");

    data.push_back(line1);
    data.push_back(line2);

    TEST_ASSERT(csv_write(data, "csvfile.csv", ";"));

    data.clear();

    TEST_ASSERT(csv_read(data, "csvfile.csv", ";"));
    TEST_EQ(2, (int)data.size());
    TEST_EQ(3, (int)data[0].size());
    TEST_EQ(4, (int)data[1].size());
    TEST_EQ(std::string("The"), data[0][0]);
    TEST_EQ(std::string("first"), data[0][1]);
    TEST_EQ(std::string("line"), data[0][2]);
    TEST_EQ(std::string("And"), data[1][0]);
    TEST_EQ(std::string("the"), data[1][1]);
    TEST_EQ(std::string("second"), data[1][2]);
    TEST_EQ(std::string("line"), data[1][3]);
    }

  }

void run_all_file_utils_tests()
  {
  using namespace jtk;
  csv_write_and_read();
  }