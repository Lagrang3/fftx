#include "wrappers.hpp"
#include <benchmark/benchmark.h>

BENCHMARK(bench_BruteForce)
    ->Arg(109)
    ->Arg(211)
    ->Arg(401)
    ->Arg(809)
    ->Arg(1009)
    ->Arg(10009)
    ->Complexity(benchmark::oNSquared);

BENCHMARK(bench_DivideAndConquer)
    ->Arg(109)
    ->Arg(211)
    ->Arg(401)
    ->Arg(809)
    ->Arg(1009)
    ->Arg(10009)
    ->Complexity(benchmark::oNSquared);

BENCHMARK(bench_Iterative)
    ->Arg(109)
    ->Arg(211)
    ->Arg(401)
    ->Arg(809)
    ->Arg(1009)
    ->Arg(10009)
    ->Complexity(benchmark::oNSquared);

BENCHMARK(bench_FFTW)
    ->Arg(109)
    ->Arg(211)
    ->Arg(401)
    ->Arg(809)
    ->Arg(1009)
    ->Arg(10009)
    ->Complexity(benchmark::oNLogN);

BENCHMARK_MAIN();
