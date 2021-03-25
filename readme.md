FFTX
===

A template library for Fast Fourier Transform.

Dependencies
===
- [Boost](https://www.boost.org/)
- [FFTW](http://fftw.org/)
- [google/benchmark](https://github.com/google/benchmark)
- [alglib](https://www.alglib.net/)
- [meson](https://mesonbuild.com/)

How to use this repository
===

This is a header only library, to use it is enough to pass a 
`-I<path to include directory>` compiler flag to your application build instructions.

The full functionality of the repository is achieved by using `meson`
buildsystem tools. For instance:
```
cd <path to build directory>
meson <path to source directory> --prefix=<path to install directory>
```
Unit tests can be run by issuing
```
meson test
```
and the library will be installed in the `<path to install directory>` once you
type
```
meson install
```

Benchmark executables will be built under `benchmarks` directory.


To-do
===
- consider the case when the algebra is non-abelian
- FFT algorithms do not depend on the value of the unity,
thus `const T _1 = T{1}` must be removed.

- Remove FFTW dependencies from this repo.
- add everything to a namespace
- use math from boost library
    https://www.boost.org/doc/libs/1_74_0/boost/algorithm/algorithm.hpp

- test the autoconf tools in an isolated environment

- parallel execution with std::thread 
    use policies to decide if run in serial or parallel
 
- benchmark against FFTW with floating point
- produce example tutorials with:
    - builtin floats
    - number theoretical fft
    - boost::multiprecision
    - boost::quaternions
- write a readme and a license
- document the use, motivation, and theory

- write primitive FFT so that we can compose FFT with those

- Unit test with boost::tests
