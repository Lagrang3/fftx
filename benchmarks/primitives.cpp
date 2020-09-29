#include <complex>
#include <functional>
#include <iostream>
#include <random>

#ifdef WITH_FFTW3
#    include <fftw3.h>
#endif

#include <benchmark/benchmark.h>

#include <fftx-1d.hpp>
#include <fftx-primitives.hpp>

typedef std::complex<double> cd;
const double PI = acos(-1.0);

std::default_random_engine gen;
std::uniform_real_distribution<double> distribution;
std::array<cd, 1000> out;

auto random_vec(size_t N)
{
    std::vector<cd> V(N);
    for (auto& x : V)
        x = distribution(gen);
    return V;
}

template <std::size_t n>
void bench_Power2(benchmark::State& state)
{
    const auto data = random_vec(n);
    for (auto _ : state)
    {
        fftx::FFT_Power2_fixed<n>(
            data.begin(), out.begin(),
            cd(cos(2 * PI / data.size()), -sin(2 * PI / data.size())));
    }
}

template <std::size_t n>
void bench_handwritten(benchmark::State& state)
{
    const auto in = random_vec(n);
    for (auto _ : state)
    {
        fftx::FFT_Handwritten_fixed<n>(in.begin(), out.begin(),
                                       cd(cos(2 * PI / n), -sin(2 * PI / n)));
    }
}
template <std::size_t n>
void bench_BruteForce(benchmark::State& state)
{
    auto data = random_vec(n);
    for (auto _ : state)
    {
        fftx::FFT_BruteForce_fixed<n>(
            data.begin(), out.begin(),
            cd(cos(2 * PI / data.size()), -sin(2 * PI / data.size())));
    }
}

#ifdef WITH_FFTW3
template <std::size_t n>
void bench_FFTW(benchmark::State& state)
{
    auto data = random_vec(n);
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * data.size());
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * data.size());

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }
    p = fftw_plan_dft_1d(data.size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    for (auto _ : state)
    {
        fftw_execute(p);
    }
    fftw_free(in);
    fftw_free(out);
    fftw_destroy_plan(p);
}
#endif

BENCHMARK_TEMPLATE(bench_Power2, 2);
BENCHMARK_TEMPLATE(bench_Power2, 4);
BENCHMARK_TEMPLATE(bench_Power2, 8);
BENCHMARK_TEMPLATE(bench_Power2, 16);
BENCHMARK_TEMPLATE(bench_Power2, 32);

BENCHMARK_TEMPLATE(bench_handwritten, 2);
BENCHMARK_TEMPLATE(bench_handwritten, 3);
BENCHMARK_TEMPLATE(bench_handwritten, 4);
BENCHMARK_TEMPLATE(bench_handwritten, 5);

BENCHMARK_TEMPLATE(bench_BruteForce, 2);
BENCHMARK_TEMPLATE(bench_BruteForce, 3);
BENCHMARK_TEMPLATE(bench_BruteForce, 4);
BENCHMARK_TEMPLATE(bench_BruteForce, 5);
BENCHMARK_TEMPLATE(bench_BruteForce, 6);
BENCHMARK_TEMPLATE(bench_BruteForce, 7);
BENCHMARK_TEMPLATE(bench_BruteForce, 8);

#ifdef WITH_FFTW3
BENCHMARK_TEMPLATE(bench_FFTW, 2);
BENCHMARK_TEMPLATE(bench_FFTW, 3);
BENCHMARK_TEMPLATE(bench_FFTW, 4);
BENCHMARK_TEMPLATE(bench_FFTW, 5);
BENCHMARK_TEMPLATE(bench_FFTW, 6);
BENCHMARK_TEMPLATE(bench_FFTW, 7);
BENCHMARK_TEMPLATE(bench_FFTW, 8);
BENCHMARK_TEMPLATE(bench_FFTW, 16);
BENCHMARK_TEMPLATE(bench_FFTW, 32);
#endif

BENCHMARK_MAIN();
