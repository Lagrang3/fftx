/*
    Demonstration of the Number Theoretical Fourier Transform.
*/

#include <algorithm>
#include <complex>
#include <fftx/1d.hpp>
#include <iostream>
#include <numeric>
#include <vector>

#include "modulo.h"

template <class T>
void print(const std::vector<T>& V)
{
    for (auto x : V)
        std::cout << x << ", ";
    std::cout << "\n";
}

int main()
{
    typedef my_modulo_lib::field_modulo<int, 337> Z337;
    using M_int = my_modulo_lib::mint<Z337>;
    const M_int w{85};
    const M_int inv_8{M_int{8}.inverse()};

    std::vector<M_int> A{4, 3, 2, 1, 0, 0, 0, 0};

    // forward FFT
    auto FT_A = fftx::FFT_Iterative(A, w);

    // backwards FFT
    auto FT_FT_A = fftx::FFT_Iterative(FT_A, w.inverse());

    std::transform(FT_FT_A.begin(), FT_FT_A.end(), FT_FT_A.begin(),
                   [&inv_8](M_int x) { return x * inv_8; });

    // compare results
    auto diff = std::transform_reduce(
        A.begin(), A.end(), FT_FT_A.begin(), 0, std::plus<M_int>(),
        [](M_int x, M_int y) { return std::abs(int(x - y)); });

    std::cout << diff << '\n';
    return 0;
}
