#include <boost/test/parameterized_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/utils/algorithm.hpp>

#include <boost/bind/bind.hpp>

#include <complex>
#include <functional>
#include <list>
#include <random>

#include <fftx-primitives.hpp>

using namespace boost::unit_test;
using namespace boost;
using namespace fftx;

// test inverse FFT

using cd = std::complex<double>;

const double PI = acos(-1.0);

template <std::size_t n, class T>
struct iterative_fft
{
    void operator()(std::vector<T>& A, const T e) const
    {
        FFT_Iterative_fixed<n>(A.begin(), e);
    }
    constexpr auto size() const { return n; }
};

template <std::size_t n, class T>
struct pow2_fft
{
    void operator()(std::vector<T>& A, const T e) const
    {
        FFT_Power2_fixed<n>(A.begin(), e);
    }
    constexpr auto size() const { return n; }
};

template <class fft_type>
void test_func_convolution(const std::vector<cd>& cA, const std::vector<cd>& cB)
{
    fft_type callable;
    BOOST_REQUIRE(cA.size() == cB.size());

    const size_t N = callable.size();
    const size_t N_half = N / 2;
    std::vector<cd> C(N, 0), A(N, 0), B(N, 0);

    std::copy(cA.begin(), cA.begin() + N_half, A.begin());
    std::copy(cB.begin(), cB.begin() + N_half, B.begin());

    // brute force convolution
    for (std::size_t i = 0; i < N; ++i)
    {
        C[i] = 0;
        for (std::size_t j = 0; j <= i; ++j)
            C[i] += A[j] * B[i - j];
    }
    // done

    // FFT-based convolution
    callable(A, cd(cos(2 * PI / N), sin(2 * PI / N)));
    callable(B, cd(cos(2 * PI / N), sin(2 * PI / N)));

    std::vector<cd> C2(N);
    for (std::size_t i = 0; i < N; ++i)
        C2[i] = A[i] * B[i];

    callable(C2, cd(cos(2 * PI / N), -sin(2 * PI / N)));

    double inv_n = 1.0 / N;
    for (auto& x : C2)
        x *= inv_n;
    // done

    // get differences
    double diff = 0;
    for (std::size_t i = 0; i < N; ++i)
    {
        diff += norm(C2[i] - C[i]);
    }
    diff = sqrt(diff) * inv_n;
    BOOST_CHECK_SMALL(diff, /*tolerance*/ 1e-8);
}

struct convolution_test_suite;

template <std::size_t beg, std::size_t end>
typename std::enable_if<beg >= end, void>::type test_linrange(
    const std::vector<cd>& A,
    const std::vector<cd>& B,
    convolution_test_suite& TS)
{
}
template <std::size_t beg, std::size_t end>
    typename std::enable_if <
    beg<end, void>::type test_linrange(const std::vector<cd>& A,
                                       const std::vector<cd>& B,
                                       convolution_test_suite& TS)
{
    TS.add(BOOST_TEST_CASE_NAME(
        std::bind(&test_func_convolution<iterative_fft<beg, cd>>, A, B),
        "Iterative fixed-size FFT, N=" + std::to_string(beg)));
    test_linrange<beg + 2, end>(A, B, TS);
}

template <std::size_t beg, std::size_t end>
typename std::enable_if<beg >= end, void>::type test_powrange(
    const std::vector<cd>& A,
    const std::vector<cd>& B,
    convolution_test_suite& TS)
{
}
template <std::size_t beg, std::size_t end>
    typename std::enable_if <
    beg<end, void>::type test_powrange(const std::vector<cd>& A,
                                       const std::vector<cd>& B,
                                       convolution_test_suite& TS)
{
    TS.add(BOOST_TEST_CASE_NAME(
        std::bind(&test_func_convolution<pow2_fft<beg, cd>>, A, B),
        "Power2 fixed-size FFT, N=" + std::to_string(beg)));
    test_powrange<beg * 2, end>(A, B, TS);
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
        // fixed size
        auto A = random_vec(128);
        auto B = random_vec(128);

        test_linrange<2, 12>(A, B, *this);
        test_powrange<2, 128>(A, B, *this);
    }
};

test_suite* init_unit_test_suite(int argc, char* argv[])
{
    framework::master_test_suite().p_name.value = "Test FFTX primitives";
    framework::master_test_suite().add(new convolution_test_suite);
    return 0;
}
