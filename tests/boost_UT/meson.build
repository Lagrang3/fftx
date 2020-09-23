test_src += [files (
    ['inverse_ut.cpp','convolution_ut.cpp'])]

boost_ut=dependency('boost',modules : ['unit_test_framework'], required: false)

if (boost_ut.found())
    convolution_ut = executable('convolution_ut',
        ['convolution_ut.cpp'],
        include_directories: [incl],
        dependencies: [boost_ut])

    test('FFT Convolution Test',convolution_ut)

    inverse_ut = executable('inverse_ut',
        ['inverse_ut.cpp'],
        include_directories: [incl],
        dependencies: [boost_ut])

    test('FFT Inverse Test',inverse_ut)
endif
