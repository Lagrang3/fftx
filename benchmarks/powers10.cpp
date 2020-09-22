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

BENCHMARK(bench_FFTW)
    ->RangeMultiplier(10)
    ->Range(10, 1000'000)
    ->Complexity(benchmark::oNLogN);

BENCHMARK_MAIN();
