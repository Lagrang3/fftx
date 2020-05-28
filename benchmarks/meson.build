test_src += [files (
    ['bench.cpp','bench_threads.cpp'])]

executable('bench',
    ['bench.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_dep])

executable('bench_threads',
    ['bench_threads.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_omp_dep,libomp,pstl,tbb])
