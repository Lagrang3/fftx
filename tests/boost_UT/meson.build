test_src += [files (
    ['inverse_ut.cpp','convolution_ut.cpp'])]


convolution_ut = executable('convolution_ut',
    ['convolution_ut.cpp'],
    include_directories: [incl],
    dependencies: [boost])

test('FFT Convolution Test',convolution_ut)

inverse_ut = executable('inverse_ut',
    ['inverse_ut.cpp'],
    include_directories: [incl],
    dependencies: [boost])

test('FFT Inverse Test',inverse_ut)
