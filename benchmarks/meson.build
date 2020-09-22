gbench = dependency('benchmark',required: false)
pthread = cpp.find_library('pthread', required: false)

fftw3 = dependency('fftw3',required: false)
fftw3_omp = dependency('fftw3_omp',required: false)
libomp=cpp.find_library('gomp', required: false)

test_src += [files (
    ['primes.cpp','powers10.cpp','powers2.cpp','wrappers.hpp'])]


if (gbench.found() and pthread.found())
    executable('powers2',['powers2.cpp'],include_directories: [incl],
        dependencies: [gbench,pthread,fftw3])
        
    executable('powers10',['powers10.cpp'],include_directories: [incl],
        dependencies: [gbench,pthread,fftw3])
        
    executable('primes',['primes.cpp'],include_directories: [incl],
        dependencies: [gbench,pthread,fftw3])
endif
