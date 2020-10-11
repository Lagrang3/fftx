#define BOOST_TEST_MODULE math
#include <boost/test/unit_test.hpp>

#include <fftx/math.hpp>

using namespace boost::unit_test;
using namespace boost;
using namespace fftx;

BOOST_AUTO_TEST_CASE(power_test)
{
    BOOST_CHECK_THROW(power(2, 0), fftx::error);

    BOOST_TEST(power(2, 1) == 2);
    BOOST_TEST(power(3, 1) == 3);

    BOOST_TEST(power(2, 2) == 4);
    BOOST_TEST(power(3, 2) == 9);

    BOOST_TEST(power(2, 3) == 8);
    BOOST_TEST(power(3, 3) == 27);

    BOOST_TEST(power(2, 4) == 16);
    BOOST_TEST(power(3, 4) == 81);

    BOOST_TEST(power(2, 5) == 32);
    BOOST_TEST(power(3, 5) == 243);

    BOOST_TEST(power(2, 6) == 64);
    BOOST_TEST(power(3, 6) == 729);
}
