#pragma once

#include <algorithm>
#include <vector>

#include <fftx-math.hpp>

/*
    Primitive Fourier transforms.

    Rules:
    - the have the signature:

    template <class iter, class T>
    void FFT_*(iter first, iter last, const T e, const T _1 = T(1));

    or

    template <std::size_t n, class iter, class T>
    void FFT_*(iter first, const T e, const T _1 = T(1));
    // here the last is not necessary because n gives the size of the container

    - they solve the FT in-place

    - iter must be a random access iterator
*/

namespace fftx
{
    /*
        handwritten fixed-size FFT
    */
    template <std::size_t n, class iter, class T>
    typename std::enable_if<n == 1, void>::type
    FFT_Handwritten_fixed(iter first, const T e, const T _1 = T(1))
    {
    }
    template <std::size_t n, class iter, class T>
    typename std::enable_if<n == 2, void>::type
    FFT_Handwritten_fixed(iter first, const T e, const T _1 = T(1))
    {
        std::array<T, n> x;
        std::copy(first, first + n, x.begin());
        first[0] = x[0] + x[1];
        first[1] = x[0] + x[1] * e;
    }
    template <std::size_t n, class iter, class T>
    typename std::enable_if<n == 3, void>::type
    FFT_Handwritten_fixed(iter first, const T e, const T _1 = T(1))
    {
        std::array<T, n> x;
        T e2 = e * e;
        std::copy(first, first + n, x.begin());
        first[0] = x[0] + x[1] + x[2];
        first[1] = x[0] + e * (x[1] + e * x[2]);
        first[2] = x[0] + e2 * (x[1] + e2 * x[2]);
    }
    template <std::size_t n, class iter, class T>
    typename std::enable_if<n == 4, void>::type
    FFT_Handwritten_fixed(iter first, const T e, const T _1 = T(1))
    {
        std::array<T, n> x, ep;
        std::copy(first, first + n, x.begin());
        // std::copy(std::execution::unsequenced_policy,first, first + n,
        // x.begin());
        ep[0] = _1;
        ep[1] = e;
        ep[2] = e * e;
        ep[3] = e * e * e;
        x[0] = first[0] + first[2];
        x[1] = first[1] + first[3];
        x[2] = first[0] + ep[2] * first[2];
        x[3] = first[1] + ep[2] * first[3];

        first[0] = x[0] + x[1];
        first[1] = x[2] + e * x[3];
        first[2] = x[0] + ep[2] * x[1];
        first[3] = x[2] + ep[3] * x[3];
    }
    template <std::size_t n, class iter, class T>
    typename std::enable_if<n == 5, void>::type
    FFT_Handwritten_fixed(iter first, const T e, const T _1 = T(1))
    {
        /*
        std::array<T, n> x;
        std::copy(first, first + n, x.begin());
        // std::copy(std::execution::unsequenced_policy,first, first + n,
        // x.begin());
        FFT_Handwritten_fixed<4>(first,e,_1);

        T e4=e*e*e*e,ep=e4;
        first[0] += x[4];
        first[1] += x[4]*ep;
        ep*=e4; first[2] += x[4]*ep;
        ep*=e4; first[3] += x[4]*ep;
        ep*=e4; first[4]=x[0]+ep*(x[1]+ep*(x[2]+ep*(x[3]+ep*x[4])));
        */
        std::array<T, n> x, ep{_1, e, e * e, e * e * e, e * e * e * e};
        std::copy(first, first + n, x.begin());
        first[0] = x[0] + x[1] + x[2] + x[3] + x[4];
        first[1] = x[0] + e * (x[1] + e * (x[2] + e * (x[3] + e * x[4])));
        first[2] =
            x[0] +
            ep[2] * (x[1] + ep[2] * (x[2] + ep[2] * (x[3] + ep[2] * x[4])));
        first[3] =
            x[0] +
            ep[3] * (x[1] + ep[3] * (x[2] + ep[3] * (x[3] + ep[3] * x[4])));
        first[4] =
            x[0] +
            ep[4] * (x[1] + ep[4] * (x[2] + ep[4] * (x[3] + ep[4] * x[4])));
    }

    /*
        fixed-size FFT with n a power of two
    */
    template <std::size_t n, class iter, class T>
    void FFT_Power2_fixed(iter first, const T e, const T _1 = T(1))
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

    /*
        fixed-size FFT with n any number
    */
    template <std::size_t n, class iter, class T>
    void FFT_Iterative_fixed(iter first, const T e, const T _1 = T(1))
    {
        std::vector<T> B(n), B_old(n);
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
