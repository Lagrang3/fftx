/*
    Simple complex to complex Discrete Fourier Transform
*/
#include <cmath>
#include <complex>
#include <fftx/1d.hpp>
#include <iostream>
#include <numeric>
#include <vector>

template <class T>
void print(const std::vector<T>& V)
{
    for (auto x : V)
        std::cout << x << ", ";
    std::cout << "\n";
}

int main()
{
    using cd = std::complex<double>;
    const cd I2Pi = cd{0.0, 1.0} * 2.0 * std::acos(-1.0);

    std::vector<cd> A{1.0, 2.0, 3.0, 4.0};
    const cd w = std::exp(I2Pi * (-1.0 / A.size()));
    print(A);

    // perform DFT on A
    auto FT_A = fftx::FFT_Iterative(A, w);
    print(FT_A);

    // DFT on FT_A with inverse root of unity (aka Inverse DFT)
    auto FT_FT_A = fftx::FFT_Iterative(FT_A, 1.0 / w);
    print(FT_FT_A);

    auto diff = std::transform_reduce(
        A.begin(), A.end(), FT_FT_A.begin(), 0.0, std::plus<double>(),
        [n = A.size()](cd x, cd y) { return std::abs(x - y * (1.0 / n)); });
    std::cout << diff << '\n';
    return 0;
}
