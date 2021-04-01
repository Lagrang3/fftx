/*
    Demonstration of the Number Theoretical Fourier Transform.
    The Schonhage-Strassen algorithm example from wikipedia
    https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm
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

    // Multiplying 1234 times 5678 = 7006652
    std::vector<M_int> A{4, 3, 2, 1, 0, 0, 0, 0};
    std::vector<M_int> B{8, 7, 6, 5, 0, 0, 0, 0};

    // forward FFT
    auto FT_A = fftx::FFT_Iterative(A, w);
    auto FT_B = fftx::FFT_Iterative(B, w);

    // convolution in Fourier space
    std::vector<M_int> FT_AB;
    std::transform(FT_A.begin(), FT_A.end(), FT_B.begin(),
                   std::back_inserter(FT_AB),
                   [](M_int x, M_int y) { return x * y; });

    // backwards FFT
    auto AB = fftx::FFT_Iterative(FT_AB, w.inverse());
    std::transform(AB.begin(), AB.end(), AB.begin(),
                   [&inv_8](M_int x) { return x * inv_8; });

    // carry the remainders in base 10
    std::vector<M_int> C;
    M_int r{0};
    for (auto x : AB)
    {
        auto y = x + r;
        C.emplace_back(int(y) % 10);
        r = M_int(int(y) / 10);
    }
    // yields 7006652
    print(C);
    return 0;
}
