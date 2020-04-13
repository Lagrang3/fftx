
fft_2d=executable('fft-2d',
    ['fft-2d.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_dep])

fft_3d=executable('fft-3d',
    ['fft-3d.cpp'],
    include_directories: [incl],
    dependencies: [fftw3_dep])

test_src += [files (
    ['fft-2d.cpp','fft-3d.cpp'])]

test('FFT-2D',fft_2d)
test('FFT-3D',fft_3d)
