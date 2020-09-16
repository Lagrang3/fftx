#include <boost/test/parameterized_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/utils/algorithm.hpp>

#include <boost/bind/bind.hpp>

#include <complex>
#include <functional>
#include <list>
#include <random>

#include <fftx.hpp>

using namespace boost::unit_test;
using namespace boost;
using namespace fftx;

// test inverse FFT

using cd = std::complex<double>;
using fft_type = decltype(FFT_BruteForce<cd>);

const double PI = acos(-1.0);

void test_func_inverse(fft_type callable, const std::vector<cd>& data)
{
    const size_t N = data.size();
    BOOST_REQUIRE(data.size() > 0);

    auto FT_data = callable(data, cd(cos(2 * PI / N), sin(2 * PI / N)), 1);
    auto FT_FT_data =
        callable(FT_data, cd(cos(2 * PI / N), -sin(2 * PI / N)), 1);

    BOOST_REQUIRE(data.size() == FT_data.size() and
                  data.size() == FT_FT_data.size());

    double inv_N = 1.0 / data.size();
    for (auto& x : FT_FT_data)
        x *= inv_N;

    double diff = 0;
    for (size_t i = 0; i < data.size(); ++i)
    {
        diff += std::norm(data[i] - FT_FT_data[i]);
    }
    diff = sqrt(diff) * inv_N;
    BOOST_CHECK_SMALL(diff, /*tolerance*/ 1e-10);
}

struct inverse_test_suite : public test_suite
{
    std::default_random_engine gen;
    std::uniform_real_distribution<double> distribution;

    auto random_vec(size_t N)
    {
        std::vector<cd> V(N);
        for (auto& x : V)
            x = distribution(gen);
        return V;
    }
    inverse_test_suite()
        : test_suite("Correctedness of the FFT, via inverse"),
          gen(123),
          distribution(0.0, 1.0)
    {
        for (size_t len = 1; len <= (1 << 12); len *= 2)
        {
            auto A = random_vec(len);

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_BruteForce<cd>, A),
                "Brute force FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_DivideAndConquer<cd>, A),
                "Divide & Conquer FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_InPlace<cd>, A),
                "In place FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_Iterative<cd>, A),
                "Iterative FFT, N=" + std::to_string(len)));
        }
        for (size_t len : {10, 100, 1000})
        {
            auto A = random_vec(len);

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_BruteForce<cd>, A),
                "Brute force FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_DivideAndConquer<cd>, A),
                "Divide & Conquer FFT, N=" + std::to_string(len)));

            // only works for powers of two
            // add(BOOST_TEST_CASE_NAME(std::bind(&test_func_inverse,
            // FFT_InPlace<cd>, A),
            //                         "In place FFT, N=" +
            //                         std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_Iterative<cd>, A),
                "Iterative FFT, N=" + std::to_string(len)));
        }
        for (size_t len : {3, 5, 7, 13, 1009})
        {
            auto A = random_vec(len);

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_BruteForce<cd>, A),
                "Brute force FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_DivideAndConquer<cd>, A),
                "Divide & Conquer FFT, N=" + std::to_string(len)));

            // only works for powers of two
            // add(BOOST_TEST_CASE_NAME(std::bind(&test_func_inverse,
            // FFT_InPlace<cd>, A),
            //                         "In place FFT, N=" +
            //                         std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_inverse, FFT_Iterative<cd>, A),
                "Iterative FFT, N=" + std::to_string(len)));
        }
    }
};

test_suite* init_unit_test_suite(int argc, char* argv[])
{
    framework::master_test_suite().p_name.value = "Boost Test FFTX";
    framework::master_test_suite().add(new inverse_test_suite);
    return 0;
}
