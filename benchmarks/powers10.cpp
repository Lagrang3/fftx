#include "wrappers.hpp"
#include <benchmark/benchmark.h>

BENCHMARK(bench_BruteForce)
    ->RangeMultiplier(10)
    ->Range(10, 10'000)
    ->Complexity(benchmark::oNSquared);

BENCHMARK(bench_DivideAndConquer)
    ->RangeMultiplier(10)
    ->Range(10, 1000'000)
    ->Complexity(benchmark::oNLogN);

BENCHMARK(bench_Iterative)
    ->RangeMultiplier(10)
    ->Range(10, 1000'000)
    ->Complexity(benchmark::oNLogN);

#ifdef WITH_FFTW3
BENCHMARK(bench_FFTW)
    ->RangeMultiplier(10)
    ->Range(10, 1000'000)
    ->Complexity(benchmark::oNLogN);
#endif

#ifdef WITH_FFTW3_OMP
BENCHMARK(bench_FFTW_omp)
    ->RangeMultiplier(10)
    ->Range(10, 1000'000)
    ->Complexity(benchmark::oNLogN);
#endif

#ifdef WITH_ALGLIB
BENCHMARK(bench_ALGLIB)
    ->RangeMultiplier(10)
    ->Range(10, 1000'000)
    ->Complexity(benchmark::oNLogN);
#endif

BENCHMARK_MAIN();
