project('fftx','cpp',
        default_options : ['cpp_std=c++17',
        'warning_level=3','optimization=3'],
        version: '0.0.1')
add_global_arguments(['-Ofast','-march=native' ],language : 'cpp')

fftw3_dep = dependency('fftw3')
fftw3_omp_dep = dependency('fftw3_omp')
cpp=meson.get_compiler('cpp')
libomp=cpp.find_library('gomp')
pstl=dependency('pstl')
tbb=dependency('tbb')

headers=[]
test_src=[]

subdir('include')
subdir('tests')
subdir('benchmarks')

_format = find_program('clang-format', required:false)
pkg=import('pkgconfig')


fftx_lib = library('fftx',headers,
    include_directories: incl,
    install: true)

install_headers(headers)

pkg.generate(fftx_lib,
    description: 'Library for FFT')

if (_format.found())
  run_target('format', command: [_format, '-i', headers,test_src])
endif
