test_src += [files (
    ['inverse_ut.cpp','convolution_ut.cpp'])]

if (boost_ut.found())
    convolution_ut = executable('convolution_ut',
        ['convolution_ut.cpp'],
        include_directories: [incl],
        dependencies: [boost_ut])

    test('FFT-primitives Convolution',convolution_ut)

    inverse_ut = executable('inverse_ut',
        ['inverse_ut.cpp'],
        include_directories: [incl],
        dependencies: [boost_ut])

    test('FFT-primitives Inverse',inverse_ut)
endif
