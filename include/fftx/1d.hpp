#pragma once

#include <algorithm>
#include <complex>
#include <fftw3.h>
#include <string>
#include <vector>

#include <fftx/math.hpp>
#include <fftx/primitives.hpp>

namespace fftx
{
    /*
        Brute force implementation of the Discrete Fourier transform.
    */
    template <class T>
    std::vector<T> FFT_BruteForce(const std::vector<T>& A,
                                  const T e,
                                  const T _1 = T(1))
    {
        const int n = A.size();
        std::vector<T> B(n);
        T ei = _1;
        for (int i = 0; i < n; ++i)
        {
            T b = 0;
            for (int j = n - 1; j >= 0; --j)
            {
                b = b * ei + A[j];
            }
            B[i] = b;
            ei *= e;
        }

        return B;
    }
    /*
        In-place FFT
        !!! n must be a power of 2 and e must be and n-root of unity (_1)
    */

    template <class iter, class T>
    void FFT_InPlace(iter first, iter last, const T e, const T _1 = T(1))
    {
        const int n = std::distance(first, last);
        if (__builtin_popcount(n) != 1)
            throw std::runtime_error(std::string(__func__) +
                                     " n=" + std::to_string(n) +
                                     " must be a power of 2");

        T f = power(e, n / 2, _1);
        int nbits = 0;
        std::vector<T> e2{e};
        for (int m = n / 2; m > 0; m >>= 1, ++nbits)
            e2.push_back(e2.back() * e2.back());

        std::reverse(e2.begin(), e2.end());

        iter _i = first;
        for (int i = 0; i < n; ++i, ++_i)
        {
            int ib = i, j = 0;
            iter _j = first;

            for (int b = 0; b < nbits; ib >>= 1, ++b)
                j = (j << 1) | (ib & 1);

            std::advance(_j, j);

            if (i < j)
                std::swap(*_i, *_j);
        }
        for (int len = 2, k = 1; len <= n; len <<= 1, ++k)
        {
            for (int i = 0; i < n; i += len)
            {
                T ej = _1;
                for (int j = 0; j < len / 2; ++j)
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
        Divide and Conquer algorithm to compute the Discrete Fourier Transform:
        aka Fast Fourier Transform.
        In a for loop
    */
    template <class T>
    std::vector<T> FFT_Iterative(const std::vector<T>& A,
                                 const T e,
                                 const T _1 = T(1))
    {
        const int n = A.size();
        std::vector<T> B(n), B_old(n);
        auto P = prime_factorization(n);

        /* reorder input  */
        for (int i = 0; i < n; ++i)
        {
            int j = 0, k = i;
            for (auto p : P)
            {
                j = j * p + k % p;
                k /= p;
            }
            B[j] = A[i];
        }

        std::reverse(P.begin(), P.end());

        /* fft */
        int len = 1;
        for (auto p : P)
        {
            int len_old = len;
            len *= p;
            std::swap(B, B_old);
            T e2 = power(e, n / len);

            for (int i = 0; i < n; i += len)
            {
                T ej = _1;
                for (int j = 0; j < len; ++j, ej *= e2)
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

        return B;
    }
    /*
        Divide and Conquer algorithm to compute the Discrete Fourier Transform:
        aka Fast Fourier Transform.
        This implementation uses recursion.
    */
    template <class T>
    std::vector<T> FFT_DivideAndConquer(const std::vector<T>& A,
                                        const T e,
                                        const T _1 = T(1))
    {
        const int n = A.size();
        if (n == 1)
            return A;

        const int p = prime_factor(n);
        const int m = n / p;
        const auto ep = power(e, p);

        std::vector<T> B(n);
        std::vector<std::vector<T>> A_sub(p, std::vector<T>(m));

        for (int i = 0; i < p; ++i)
        {
            for (int j = 0; j < m; ++j)
                A_sub[i][j] = A[j * p + i];

            A_sub[i] = FFT_DivideAndConquer(A_sub[i], ep, _1);
        }

        T ek = _1;
        for (int k = 0; k < n; ++k)
        {
            T b = 0;
            for (int i = p - 1; i >= 0; --i)
            {
                b = b * ek + A_sub[i][k % m];
            }
            B[k] = b;
            ek *= e;
        }

        return B;
    }

    /*
        Divide and Conquer algorithm to compute the Discrete Fourier Transform:
        aka Fast Fourier Transform.
        In a for loop

        TBB threaded version
    */
    /*
    template <class T>
    std::vector<T> FFT_Iterative_parallel(const std::vector<T>& A,
                                          const T e,
                                          const T _1 = T(1))
    {
        const int n = A.size();
        std::vector<T> B(n), B_old(n);
        auto P = prime_factorization(n);

        // reorder input
        tbb::parallel_for(0, n, 100, [&A, &B, &P](int i) {
            int j = 0, k = i;
            for (auto p : P)
            {
                j = j * p + k % p;
                k /= p;
            }
            B[j] = A[i];
        });

        std::reverse(P.begin(), P.end());

        // fft
        int len = 1;
        for (auto p : P)
        {
            int len_old = len;
            len *= p;
            std::swap(B, B_old);
            T e2 = power(e, n / len);

            tbb::parallel_for(0, n / len, 1, [&](int idx) {
                int i = idx * len;
                T ej = _1;
                for (int j = 0; j < len; ++j, ej *= e2)
                {
                    if (p < 2000)
                    {
                        T b = 0;
                        for (int k = p - 1; k >= 0; --k)
                        {
                            b = b * ej + B_old[i + k * len_old + j % len_old];
                        }
                        B[i + j] = b;
                    }
                    else
                    {
                        using pair = std::pair<T, int>;

                        B[i + j] =
                            tbb::parallel_reduce(
                                tbb::blocked_range<int>(0, p, 100),
                                pair(T(0), 0),
                                [&](const tbb::blocked_range<int>& range,
                                    pair init) -> pair {
                                    T b = 0;
                                    for (int k = range.end() - 1;
                                         k >= range.begin(); --k)
                                    {
                                        b = b * ej + B_old[i + k * len_old +
                                                           j % len_old];
                                    }
                                    return pair(b, range.size());
                                },
                                [&ej](pair a, pair b) {
                                    return pair(
                                        a.first + b.first * power(ej, a.second),
                                        a.second + b.second);
                                })
                                .first;
                    }
                }
            });
        }

        return B;
    }
    */
    /*
        Divide and Conquer algorithm to compute the Discrete Fourier Transform:
        aka Fast Fourier Transform.
        !!! n must be a power of 2 and e must be and n-root of unity (_1)

        In Place wrapper
    */
    template <class T>
    std::vector<T> FFT_InPlace(const std::vector<T>& A,
                               const T e,
                               const T _1 = T(1))
    {
        std::vector<T> B(A);
        FFT_InPlace(B.begin(), B.end(), e, _1);
        return B;
    }

}  // namespace fftx
