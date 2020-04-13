project('fftx','cpp',
        default_options : ['cpp_std=c++17',
        'warning_level=3','optimization=3'],
        version: '0.0.1')

fftw3_dep = dependency('fftw3')

headers=[]
test_src=[]

subdir('include')
subdir('tests')

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
