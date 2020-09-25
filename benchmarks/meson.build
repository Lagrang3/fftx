gbench = dependency('benchmark',required: false)
pthread = cpp.find_library('pthread', required: false)

fftw3 = dependency('fftw3',required: false)
fftw3_omp = dependency('fftw3_omp',required: false)
omp=cpp.find_library('gomp', required: false)

alglib=dependency('alglib', required: false)

test_src += [files (
    ['primes.cpp','powers10.cpp','powers2.cpp','wrappers.hpp',
    'primitives.cpp'])]


if (gbench.found() and pthread.found())
    dep = [gbench,pthread]
    args = []
    include = [incl,'.']
    
    if (fftw3.found())
        dep += [fftw3]
        args += ['-DWITH_FFTW3']
    endif
    
    if (fftw3_omp.found() and omp.found())
        dep += [omp,fftw3_omp]
        args += ['-DWITH_FFTW3_OMP']
    endif
    
    if (alglib.found())
        dep += [alglib]
        args += ['-DWITH_ALGLIB']
    endif
    
    executable('powers2',['powers2.cpp'],include_directories: include,
        dependencies: dep,
        cpp_args: args)
        
    executable('powers10',['powers10.cpp'],include_directories: include,
        dependencies: dep,
        cpp_args: args)
        
    executable('primes',['primes.cpp'],include_directories: include,
        dependencies: dep,
        cpp_args: args)
    
    executable('primitives',['primitives.cpp'],include_directories: include,
        dependencies: dep,
        cpp_args: args)
endif
