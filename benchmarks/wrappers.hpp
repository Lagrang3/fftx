#include <complex>
#include <functional>
#include <random>

#ifdef WITH_FFTW3
#    include <fftw3.h>
#endif

#include <benchmark/benchmark.h>

#include <fftx.hpp>

typedef std::complex<double> cd;
const double PI = acos(-1.0);
using fft_type = decltype(fftx::FFT_BruteForce<cd>);

std::default_random_engine gen;
std::uniform_real_distribution<double> distribution;

auto random_vec(size_t N)
{
    std::vector<cd> V(N);
    for (auto& x : V)
        x = distribution(gen);
    return V;
}

void bench_BruteForce(benchmark::State& state)
{
    auto data = random_vec(state.range(0));
    for (auto _ : state)
    {
        fftx::FFT_BruteForce<cd>(
            data, cd(cos(2 * PI / data.size()), -sin(2 * PI / data.size())), 1);
    }
    state.SetComplexityN(state.range(0));
}

void bench_InPlace(benchmark::State& state)
{
    auto data = random_vec(state.range(0));
    for (auto _ : state)
    {
        fftx::FFT_InPlace<cd>(
            data, cd(cos(2 * PI / data.size()), -sin(2 * PI / data.size())), 1);
    }
    state.SetComplexityN(state.range(0));
}

void bench_DivideAndConquer(benchmark::State& state)
{
    auto data = random_vec(state.range(0));
    for (auto _ : state)
    {
        fftx::FFT_DivideAndConquer<cd>(
            data, cd(cos(2 * PI / data.size()), -sin(2 * PI / data.size())), 1);
    }
    state.SetComplexityN(state.range(0));
}
void bench_Iterative(benchmark::State& state)
{
    auto data = random_vec(state.range(0));
    for (auto _ : state)
    {
        fftx::FFT_Iterative<cd>(
            data, cd(cos(2 * PI / data.size()), -sin(2 * PI / data.size())), 1);
    }
    state.SetComplexityN(state.range(0));
}

#ifdef WITH_FFTW3
void bench_FFTW(benchmark::State& state)
{
    auto data = random_vec(state.range(0));
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * data.size());
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * data.size());

    for (size_t i = 0; i < data.size(); ++i)
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
    state.SetComplexityN(state.range(0));
}
#endif

#ifdef WITH_FFTW3_OMP
void bench_FFTW_omp(benchmark::State& state)
{
    fftw_init_threads();
    fftw_plan_with_nthreads(12);
    auto data = random_vec(state.range(0));
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * data.size());
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * data.size());

    for (size_t i = 0; i < data.size(); ++i)
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
    state.SetComplexityN(state.range(0));
    fftw_cleanup_threads();
}
#endif
