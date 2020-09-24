#include "wrappers.hpp"
#include <benchmark/benchmark.h>

BENCHMARK(bench_BruteForce)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 15)
    ->Complexity(benchmark::oNSquared);

BENCHMARK(bench_InPlace)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

BENCHMARK(bench_DivideAndConquer)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

BENCHMARK(bench_Iterative)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);

#ifdef WITH_FFTW3
BENCHMARK(bench_FFTW)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);
#endif

#ifdef WITH_FFTW3_OMP
BENCHMARK(bench_FFTW_omp)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);
#endif

#ifdef WITH_ALGLIB
BENCHMARK(bench_ALGLIB)
    ->RangeMultiplier(8)
    ->Range(1 << 5, 1 << 20)
    ->Complexity(benchmark::oNLogN);
#endif

BENCHMARK_MAIN();
