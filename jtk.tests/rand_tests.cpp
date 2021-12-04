#include "rand_tests.h"
#include "test_assert.h"
#include "../jtk/rand.h"

#include <iostream>
#include <random>

namespace jtk
  {


  void test_xorshift32()
    {
    xorshift32 gen;

    gen.seed(0x74382381);

    TEST_EQ(2561652095, gen());
    TEST_EQ(2920033695, gen());
    TEST_EQ(1548013632, gen());
    TEST_EQ(266002342, gen());
    TEST_EQ(3521676177, gen());
    TEST_EQ(4072843484, gen());
    TEST_EQ(1559633200, gen());
    TEST_EQ(1968832920, gen());
    TEST_EQ(2555018924, gen());
    TEST_EQ(611137795, gen());    
    }

  void test_xorshift64()
    {
    xorshift64 gen;

    gen.seed(0x74382381);

    TEST_EQ(7967805708226634297, gen());
    TEST_EQ(9873696482314704701, gen());
    TEST_EQ(1707834547925741779, gen());
    TEST_EQ(3432389062673704434, gen());
    TEST_EQ(277581967194074425, gen());
    TEST_EQ(15535706225997976355, gen());
    TEST_EQ(1474474744314815941, gen());
    TEST_EQ(5716425719980810830, gen());
    TEST_EQ(11805380213173979146, gen());
    TEST_EQ(1255979623467391130, gen());
    }

  void test_xorshift64star()
    {
    xorshift64star gen;

    gen.seed(0x74382381);

    TEST_EQ(18002989874742077248, gen());
    TEST_EQ(16746940525103911697, gen());
    TEST_EQ(6262314654144044813, gen());
    TEST_EQ(11174128747390327810, gen());
    TEST_EQ(14090191277248588175, gen());
    TEST_EQ(18229059219140614279, gen());
    TEST_EQ(16392023110057125236, gen());
    TEST_EQ(1440475888802347666, gen());
    TEST_EQ(11424715604177284307, gen());
    TEST_EQ(9695351890346187618, gen());
    }

  void test_rndhash32()
    {
    rndhash32 gen;    

    TEST_EQ(436901570, gen());
    TEST_EQ(725778245, gen());
    TEST_EQ(61355516, gen());
    TEST_EQ(3809932532, gen());
    TEST_EQ(2923951477, gen());
    TEST_EQ(4286858955, gen());
    TEST_EQ(357163824, gen());
    TEST_EQ(2278085025, gen());
    TEST_EQ(1639418064, gen());
    TEST_EQ(2704673657, gen());
    }

  void test_rndhash64()
    {
    rndhash64 gen;

    TEST_EQ(3552937547537979953, gen());
    TEST_EQ(6794055111782360197, gen());
    TEST_EQ(5856732661605880733, gen());
    TEST_EQ(10810044539420750887, gen());
    TEST_EQ(15310489133201607694, gen());
    TEST_EQ(10529180171964038018, gen());
    TEST_EQ(5631119354899412717, gen());
    TEST_EQ(432025260575737302, gen());
    TEST_EQ(12448303717905564719, gen());
    TEST_EQ(9360494376432242552, gen());
    }

  void test_uniform_int_distribution()
    {
    rndhash32 gen;
    std::uniform_int_distribution<int> dist(0,100);
    TEST_EQ(12, dist(gen));
    xorshift32 gen2;
    TEST_EQ(11, dist(gen2));
    xorshift64 gen3;
    TEST_EQ(22, dist(gen3));
    xorshift64star gen4;
    TEST_EQ(33, dist(gen4));
    rndhash64 gen5;
    TEST_EQ(95, dist(gen5));
    }

  void test_uniform_real_distribution()
    {
    rndhash32 gen;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    TEST_EQ(0.16898341595701036, dist(gen));
    xorshift32 gen2;
    TEST_EQ(0.52092439366361742, dist(gen2));
    xorshift64 gen3;
    TEST_EQ(0.53525415882940430, dist(gen3));
    xorshift64star gen4;
    TEST_EQ(0.69168813734950074, dist(gen4));
    rndhash64 gen5;
    TEST_EQ(0.19260513038729990, dist(gen5));
    }
  }

void run_all_rand_tests()
  {
  using namespace jtk;
  test_xorshift32();
  test_xorshift64();
  test_xorshift64star();
  test_rndhash32();
  test_rndhash64();
  test_uniform_int_distribution();
  test_uniform_real_distribution();
  }
