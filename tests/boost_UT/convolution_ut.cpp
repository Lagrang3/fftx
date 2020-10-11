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

using cd = std::complex<double>;
using fft_type = decltype(FFT_BruteForce<cd>);

const double PI = acos(-1.0);

// test with convolution
void test_func_convolution(fft_type callable,
                           const std::vector<cd>& cA,
                           const std::vector<cd>& cB)
{
    BOOST_REQUIRE(cA.size() > 0);
    BOOST_REQUIRE(cA.size() == cB.size());

    const int N = cA.size() * 2;
    std::vector<cd> C(N), A{cA}, B{cB};

    A.resize(N, 0);
    B.resize(N, 0);

    // brute force convolution
    for (int i = 0; i < N; ++i)
    {
        C[i] = 0;
        for (int j = 0; j <= i; ++j)
            C[i] += A[j] * B[i - j];
    }
    // done

    // FFT-based convolution
    auto FA = callable(A, cd(cos(2 * PI / N), sin(2 * PI / N)));
    auto FB = callable(B, cd(cos(2 * PI / N), sin(2 * PI / N)));

    std::vector<cd> FC(N);
    for (int i = 0; i < N; ++i)
        FC[i] = FA[i] * FB[i];

    auto C2 = callable(FC, cd(cos(2 * PI / N), -sin(2 * PI / N)));

    double inv_n = 1.0 / N;
    for (auto& x : C2)
        x *= inv_n;
    // done

    // get differences
    double diff = 0;
    for (int i = 0; i < N; ++i)
    {
        diff += norm(C2[i] - C[i]);
    }
    diff = sqrt(diff) * inv_n;
    BOOST_CHECK_SMALL(diff, /*tolerance*/ 1e-8);
}
struct convolution_test_suite : public test_suite
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
    convolution_test_suite()
        : test_suite("Correctedness of the FFT, via convolution"),
          gen(123),
          distribution(0.0, 1.0)
    {
        for (size_t len = 1; len <= (1 << 12); len *= 2)
        {
            auto A = random_vec(len);
            auto B = random_vec(len);

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_BruteForce<cd>, A, B),
                "Brute force FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_DivideAndConquer<cd>, A,
                          B),
                "Divide & Conquer FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_InPlace<cd>, A, B),
                "In place FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_Iterative<cd>, A, B),
                "Iterative FFT, N=" + std::to_string(len)));
        }
        for (size_t len : {10, 100, 1000})
        {
            auto A = random_vec(len);
            auto B = random_vec(len);

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_BruteForce<cd>, A, B),
                "Brute force FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_DivideAndConquer<cd>, A,
                          B),
                "Divide & Conquer FFT, N=" + std::to_string(len)));

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_Iterative<cd>, A, B),
                "Iterative FFT, N=" + std::to_string(len)));
        }
        for (size_t len : {3, 5, 7, 13, 1009})
        {
            auto A = random_vec(len);
            auto B = random_vec(len);

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_BruteForce<cd>, A, B),
                "Brute force FFT, N=" + std::to_string(len)));
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_DivideAndConquer<cd>, A,
                          B),
                "Divide & Conquer FFT, N=" + std::to_string(len)));

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_func_convolution, FFT_Iterative<cd>, A, B),
                "Iterative FFT, N=" + std::to_string(len)));
        }
    }
};

test_suite* init_unit_test_suite(int argc, char* argv[])
{
    framework::master_test_suite().p_name.value = "Boost Test FFTX";
    framework::master_test_suite().add(new convolution_test_suite);
    return 0;
}
