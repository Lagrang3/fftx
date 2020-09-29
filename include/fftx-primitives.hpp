#pragma once

#include <algorithm>
#include <vector>

#include <fftx-math.hpp>

/*
    Primitive Fourier transforms.

    Rules:
    - the have the signature:

    template <class iter, class T>
    void FFT_*(iter first, iter last, const T e);

    or

    template <std::size_t n, class iter, class T>
    void FFT_*(iter first, const T e);
    // here the last is not necessary because n gives the size of the container

    - they solve the FT in-place

    - iter must be a random access iterator
*/

namespace fftx
{
    /*
        handwritten fixed-size FFT
    */
    template <std::size_t n, class iter1, class iter2, class T>
    typename std::enable_if<n == 1, void>::type FFT_Handwritten_fixed(iter1 in,
                                                                      iter2 out,
                                                                      const T e)
    {
        out[0] = in[0];
    }

    template <std::size_t n, class iter1, class iter2, class T>
    typename std::enable_if<n == 2, void>::type FFT_Handwritten_fixed(iter1 in,
                                                                      iter2 out,
                                                                      const T e)
    {
        // we use x in case 'in' and 'out' point to the same location
        std::array<T, n> x;
        std::copy(in, in + n, x.begin());
        out[0] = x[0] + x[1];
        out[1] = x[0] + e * x[1];
    }

    template <std::size_t n, class iter1, class iter2, class T>
    typename std::enable_if<n == 3, void>::type FFT_Handwritten_fixed(iter1 in,
                                                                      iter2 out,
                                                                      const T e)
    {
        // we use x in case 'in' and 'out' point to the same location
        std::array<T, n> x;
        std::copy(in, in + n, x.begin());
        T e2 = e * e;
        out[0] = x[0] + x[1] + x[2];
        out[1] = x[0] + e * x[1] + e2 * x[2];
        out[2] = x[0] + e2 * x[1] + e2 * e2 * x[2];
    }

    template <std::size_t n, class iter1, class iter2, class T>
    typename std::enable_if<n == 4, void>::type FFT_Handwritten_fixed(iter1 in,
                                                                      iter2 out,
                                                                      const T e)
    {
        std::array<T, n> x;
        T e2 = e * e, e3 = e2 * e;
        x[0] = in[0] + in[2];
        x[1] = in[1] + in[3];
        x[2] = in[0] + e2 * in[2];
        x[3] = in[1] + e2 * in[3];

        out[0] = x[0] + x[1];
        out[1] = x[2] + e * x[3];
        out[2] = x[0] + e2 * x[1];
        out[3] = x[2] + e3 * x[3];
    }

    template <std::size_t n, class iter1, class iter2, class T>
    typename std::enable_if<n == 5, void>::type FFT_Handwritten_fixed(iter1 in,
                                                                      iter2 out,
                                                                      const T e)
    {
        std::array<T, n> x;
        T e2 = e * e, e3 = e2 * e, e4 = e3 * e;
        std::copy(in, in + n, x.begin());

        out[0] = x[0] + x[1] + x[2] + x[3] + x[4];
        for (std::size_t i = 1; i < n; ++i)
        {
            x[1] *= e;
            x[2] *= e2;
            x[3] *= e3;
            x[4] *= e4;
            out[i] = x[0] + x[1] + x[2] + x[3] + x[4];
        }
    }

    template <std::size_t n, class iter1, class iter2, class T>
    void FFT_BruteForce_fixed(iter1 in, iter2 out, const T e)
    {
        if (n == 1)
            out[0] = in[0];
        std::array<T, n> x;
        std::copy(in, in + n, x.begin());

        // we cannot assume that T{0} is the null element
        // we actually don't need the null element
        out[0] = x[0];
        for (std::size_t i = 1; i < n; ++i)
            out[0] += x[i];

        T ei = e;
        for (std::size_t i = 1; i < n; ++i)
        {
            T b = x[n - 1];
            for (int j = n - 2; j >= 0; --j)
            {
                b = b * ei + x[j];
            }
            out[i] = b;
            ei *= e;
        }
    }

    /*
        fixed-size FFT with n a power of two
    */
    template <std::size_t n, class iter1, class iter2, class T>
    void FFT_Power2_fixed(iter1 in, iter2 out, const T e)
    {
        if (n == 1)
            out[0] = in[0];
        std::array<T, n> x;
        std::copy(in, in + n, x.begin());

        T f = power(e, n / 2);
        int nbits = 0;
        std::vector<T> e2{e};
        for (int m = n / 2; m > 0; m >>= 1, ++nbits)
            e2.push_back(e2.back() * e2.back());

        std::reverse(e2.begin(), e2.end());

        for (std::size_t i = 0; i < n; ++i)
        {
            std::size_t ib = i, j = 0;

            for (int b = 0; b < nbits; ib >>= 1, ++b)
                j = (j << 1) | (ib & 1);

            if (i < j)
                std::swap(x[i], x[j]);
        }
        for (std::size_t len = 2, k = 1; len <= n; len <<= 1, ++k)
        {
            for (std::size_t i = 0; i < n; i += len)
            {
                {
                    // j=0
                    T &u = x[i], &v = x[i + len / 2];
                    T Bu = u, Bv = v;
                    u = Bu + Bv;
                    v = Bu + Bv * f;
                }
                T ej = e2[k];
                for (std::size_t j = 1; j < len / 2; ++j)
                {
                    T &u = x[i + j], &v = x[i + j + len / 2];
                    T Bu = u, Bv = v * ej;
                    u = Bu + Bv;
                    v = Bu + Bv * f;
                    ej *= e2[k];
                }
            }
        }
        std::copy(x.begin(), x.end(), out);
    }

    /*
        fixed-size FFT with n any number
    */
    template <std::size_t n, class iter1, class iter2, class T>
    void FFT_Iterative_fixed(iter1 in, iter2 out, const T e)
    {
        if (n == 1)
            out[0] = in[0];
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
            B[j] = in[i];
        }

        std::reverse(P.begin(), P.end());

        /* fft */
        std::size_t len = 1;
        for (auto p : P)
        {
            int len_old = len;
            len *= p;
            std::swap(B, B_old);
            T e2 = power(e, n / len);  // len<=n and len divides n

            for (std::size_t i = 0; i < n; i += len)
            {
                {
                    // j =0
                    T b = B_old[i + (p - 1) * len_old];
                    for (int k = p - 2; k >= 0; --k)
                    {
                        b += B_old[i + k * len_old];
                    }
                    B[i] = b;
                }
                T ej = e2;
                for (std::size_t j = 1; j < len; ++j, ej *= e2)
                {
                    T b = B_old[i + (p - 1) * len_old + j % len_old];
                    for (int k = p - 2; k >= 0; --k)
                    {
                        b = b * ej + B_old[i + k * len_old + j % len_old];
                    }
                    B[i + j] = b;
                }
            }
        }

        std::copy(B.begin(), B.end(), out);
    }
}  // namespace fftx
