test_src += [files (
    ['bench.cpp','bench_cuda.cpp','fft.cu'])]

executable('bench',
    ['bench.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_dep])

executable('bench_cuda',
    'bench_cuda.cpp','fft.cu',
    include_directories: [incl],
    dependencies: [cufft])
