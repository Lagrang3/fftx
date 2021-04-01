example_include = include_directories('.')
example_src += [files (
    ['ex01.cpp','ex02.cpp','ex03.cpp','modulo.h'])]

ex01 = executable('ex01',
   ['ex01.cpp'],
   include_directories: [incl,example_include])
ex02 = executable('ex02',
   ['ex02.cpp','modulo.h'],
   include_directories: [incl,example_include])
ex03 = executable('ex03',
   ['ex03.cpp','modulo.h'],
   include_directories: [incl,example_include])
