test_src += [files (
    ['inverse_ut.cpp','convolution_ut.cpp','math.cpp'])]

if (boost_ut.found())
   convolution_ut = executable('convolution_ut',
       ['convolution_ut.cpp'],
       include_directories: [incl],
       dependencies: [boost_ut])

   test('FFT Convolution',convolution_ut)

    inverse_ut = executable('inverse_ut',
        ['inverse_ut.cpp'],
        include_directories: [incl],
        dependencies: [boost_ut])

    test('FFT Inverse',inverse_ut)
    
    math_ut = executable('math_ut',
        ['math.cpp'],
        include_directories: [incl],
        dependencies: [boost_ut])

    test('Math',math_ut)
endif
