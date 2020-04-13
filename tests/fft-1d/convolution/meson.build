brute_force=executable('brute-force',
    ['brute-force.cpp'],
    include_directories: [incl])

divide_n_conquer=executable('divide-n-conquer',
    ['divide-n-conquer.cpp'],
    include_directories: [incl])

in_place=executable('in-place',
    ['in-place.cpp'],
    include_directories: [incl])

test_src += [files (
    ['brute-force.cpp','divide-n-conquer.cpp','in-place.cpp'])]

generator = find_program('gen_test.py')

test_files = custom_target('test-files',
  output : ['01.in','02.in','03.in','04.in','05.in'],
  command : [generator,'@OUTDIR@'])

foreach _t : test_files.to_list()
    test('FFT-1D Convolution brute-force',
        judge,
        args: [brute_force.full_path(),_t.full_path()],
        depends: test_files)
endforeach
    
foreach _t : test_files.to_list()
    test('FFT-1D Convolution divide-n-conquer',
        judge,
        args: [divide_n_conquer.full_path(),_t.full_path()],
        depends: test_files)
endforeach
    
foreach _t : test_files.to_list()
    test('FFT-1D Convolution in-place',
        judge,
        args: [in_place.full_path(),_t.full_path()],
        depends: test_files)
endforeach
