test_src += [files (
    ['bench.cpp','bench_threads.cpp',
        'primes.cpp','powers10.cpp','powers2.cpp','wrappers.hpp'])]

executable('bench',
    ['bench.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_dep])

executable('bench_threads',
    ['bench_threads.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_omp_dep,libomp,pstl,tbb])

gbench = dependency('benchmark')
pthread = cpp.find_library('pthread')

executable('powers2',['powers2.cpp'],include_directories: [incl],
    dependencies: [gbench,pthread,fftw3_dep])
    
executable('powers10',['powers10.cpp'],include_directories: [incl],
    dependencies: [gbench,pthread,fftw3_dep])
    
executable('primes',['primes.cpp'],include_directories: [incl],
    dependencies: [gbench,pthread,fftw3_dep])
