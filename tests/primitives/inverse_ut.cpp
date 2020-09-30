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
        FFT_Iterative_fixed<n>(A.begin(), A.begin(), e);
    }
    constexpr auto size() const { return n; }
};

template <std::size_t n, class T>
struct pow2_fft
{
    void operator()(std::vector<T>& A, const T e) const
    {
        FFT_Power2_fixed<n>(A.begin(), A.begin(), e);
    }
    constexpr auto size() const { return n; }
};

template <std::size_t n, class T>
struct handwritten_fft
{
    void operator()(std::vector<T>& A, const T e) const
    {
        FFT_Handwritten_fixed<n>(A.begin(), A.begin(), e);
    }
    constexpr auto size() const { return n; }
};

template <std::size_t n, class T>
struct bruteforce_fft
{
    void operator()(std::vector<T>& A, const T e) const
    {
        FFT_BruteForce_fixed<n>(A.begin(), A.begin(), e);
    }
    constexpr auto size() const { return n; }
};

template <class fft_type>
void test_func_inverse(const std::vector<cd>& data)
{
    fft_type callable;
    const size_t N = callable.size();

    std::vector<cd> B(N);

    std::copy(data.begin(), data.begin() + N, B.begin());

    callable(B, cd(cos(2 * PI / N), sin(2 * PI / N)));
    callable(B, cd(cos(2 * PI / N), -sin(2 * PI / N)));

    double inv_N = 1.0 / N;
    for (auto& x : B)
        x *= inv_N;

    double diff = 0;
    for (size_t i = 0; i < N; ++i)
    {
        diff += std::norm(data[i] - B[i]);
    }
    diff = sqrt(diff) * inv_N;
    BOOST_CHECK_SMALL(diff, /*tolerance*/ 1e-10);
}

struct inverse_test_suite;

template <std::size_t beg, std::size_t end>
typename std::enable_if<beg >= end, void>::type test_linrange(
    const std::vector<cd>& A,
    inverse_test_suite& TS)
{
}
template <std::size_t beg, std::size_t end>
    typename std::enable_if <
    beg<end, void>::type test_linrange(const std::vector<cd>& A,
                                       inverse_test_suite& TS)
{
    TS.add(BOOST_TEST_CASE_NAME(
        std::bind(&test_func_inverse<iterative_fft<beg, cd>>, A),
        "Iterative fixed-size FFT, N=" + std::to_string(beg)));
    test_linrange<beg + 1, end>(A, TS);
}

template <std::size_t beg, std::size_t end>
typename std::enable_if<beg >= end, void>::type test_powrange(
    const std::vector<cd>& A,
    inverse_test_suite& TS)
{
}
template <std::size_t beg, std::size_t end>
    typename std::enable_if <
    beg<end, void>::type test_powrange(const std::vector<cd>& A,
                                       inverse_test_suite& TS)
{
    TS.add(BOOST_TEST_CASE_NAME(
        std::bind(&test_func_inverse<pow2_fft<beg, cd>>, A),
        "Power2 fixed-size FFT, N=" + std::to_string(beg)));
    test_powrange<beg * 2, end>(A, TS);
}

template <std::size_t beg, std::size_t end>
typename std::enable_if<beg >= end, void>::type test_linrange_handwritten(
    const std::vector<cd>& A,
    inverse_test_suite& TS)
{
}
template <std::size_t beg, std::size_t end>
    typename std::enable_if <
    beg<end, void>::type test_linrange_handwritten(const std::vector<cd>& A,
                                                   inverse_test_suite& TS)
{
    TS.add(BOOST_TEST_CASE_NAME(
        std::bind(&test_func_inverse<handwritten_fft<beg, cd>>, A),
        "Handwritten fixed-size FFT, N=" + std::to_string(beg)));
    test_linrange_handwritten<beg + 1, end>(A, TS);
}

template <std::size_t beg, std::size_t end>
typename std::enable_if<beg >= end, void>::type test_linrange_bruteforce(
    const std::vector<cd>& A,
    inverse_test_suite& TS)
{
}
template <std::size_t beg, std::size_t end>
    typename std::enable_if <
    beg<end, void>::type test_linrange_bruteforce(const std::vector<cd>& A,
                                                  inverse_test_suite& TS)
{
    TS.add(BOOST_TEST_CASE_NAME(
        std::bind(&test_func_inverse<bruteforce_fft<beg, cd>>, A),
        "BruteForce fixed-size FFT, N=" + std::to_string(beg)));
    test_linrange_bruteforce<beg + 1, end>(A, TS);
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
        // fixed size
        auto A = random_vec(128);

        test_linrange<1, 10>(A, *this);
        test_powrange<1, 128>(A, *this);

        test_linrange_handwritten<1, 8>(A, *this);
        test_linrange_bruteforce<1, 10>(A, *this);
    }
};

test_suite* init_unit_test_suite(int argc, char* argv[])
{
    framework::master_test_suite().p_name.value = "Test FFTX primitives";
    framework::master_test_suite().add(new inverse_test_suite);
    return 0;
}
