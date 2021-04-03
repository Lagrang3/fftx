/*
    Demonstration of the Number Theoretical Fourier Transform.
*/

#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
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

class Z337
{
   public:
    typedef boost::multiprecision::int128_t integer;
    static constexpr integer mod = 337;
};

int main()
{
    // using namespace boost::multiprecision;

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
    int diff = 0;
    for (size_t i = 0; i < A.size(); ++i)
        diff += A[i] == FT_FT_A[i] ? 0 : 1;

    std::cout << diff << '\n';
    return 0;
}
