project('fftx','cpp','cuda',
        version: '0.0.1')
        
add_global_arguments([
    '-O3','-Ofast','-march=native','-std=c++17','-Wall'],
    language : 'cpp')
    
add_global_arguments(['-O3'],language : 'cuda')

fftw3_dep = dependency('fftw3')
cufft = dependency('cuda', modules: ['cufft'] )

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
