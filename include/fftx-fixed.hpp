#pragma once
#include <fftx-math.hpp>

/*
    Fixed size Fourier transforms
*/

namespace fftx
{
    template <std::size_t n, class iter, class T>
    void FFT_InPlace_fixed(iter first, const T e, const T _1 = T(1))
    {
        T f = power(e, n / 2, _1);
        int nbits = 0;
        std::vector<T> e2{e};
        for (int m = n / 2; m > 0; m >>= 1, ++nbits)
            e2.push_back(e2.back() * e2.back());

        std::reverse(e2.begin(), e2.end());

        iter _i = first;
        for (std::size_t i = 0; i < n; ++i, ++_i)
        {
            std::size_t ib = i, j = 0;
            iter _j = first;

            for (int b = 0; b < nbits; ib >>= 1, ++b)
                j = (j << 1) | (ib & 1);

            std::advance(_j, j);

            if (i < j)
                std::swap(*_i, *_j);
        }
        for (std::size_t len = 2, k = 1; len <= n; len <<= 1, ++k)
        {
            for (std::size_t i = 0; i < n; i += len)
            {
                T ej = _1;
                for (std::size_t j = 0; j < len / 2; ++j)
                {
                    iter u = first + i + j, v = first + i + j + len / 2;
                    T Bu = *u, Bv = *v * ej;
                    *u = Bu + Bv;
                    *v = Bu + Bv * f;
                    ej *= e2[k];
                }
            }
        }
    }

    template <std::size_t n, class iter, class T>
    void FFT_Iterative_fixed(iter first, const T e, const T _1 = T(1))
    {
        std::array<T, n> B, B_old;
        auto P = prime_factorization(n);

        /* reorder input  */
        for (std::size_t i = 0; i < n; ++i)
        {
            int j = 0, k = i;
            for (auto p : P)
            {
                j = j * p + k % p;
                k /= p;
            }
            B[j] = first[i];
        }

        std::reverse(P.begin(), P.end());

        /* fft */
        std::size_t len = 1;
        for (auto p : P)
        {
            int len_old = len;
            len *= p;
            std::swap(B, B_old);
            T e2 = power(e, n / len);

            for (std::size_t i = 0; i < n; i += len)
            {
                T ej = _1;
                for (std::size_t j = 0; j < len; ++j, ej *= e2)
                {
                    T b = 0;
                    for (int k = p - 1; k >= 0; --k)
                    {
                        b = b * ej + B_old[i + k * len_old + j % len_old];
                    }
                    B[i + j] = b;
                }
            }
        }

        std::copy(B.begin(), B.end(), first);
    }
}  // namespace fftx
