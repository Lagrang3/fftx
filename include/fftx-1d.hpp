#pragma once

#include <algorithm>
#include <cassert>
#include <complex>
#include <fftw3.h>
#include <vector>

#include <fftx-math.hpp>

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
        T eij = _1;
        B[i] = 0;
        for (int j = 0; j < n; ++j)
        {
            B[i] += A[j] * eij;
            eij *= ei;
        }
        ei *= e;
    }

    return B;
}

/*
    Divide and Conquer algorithm to compute the Discrete Fourier Transform:
    aka Fast Fourier Transform.
    This implementation uses recursion.
    !!! n must be a power of 2 and e must be and n-root of unity (_1)
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

    if (p == n)
        return FFT_BruteForce(A, e, _1);

    const int m = n / p;
    const auto ep = power(e, p);

    std::vector<T> B(n);
    std::vector<std::vector<T>> A_sub(p);

    for (int i = 0; i < p; ++i)
    {
        std::vector<T> tmp(m);
        for (int j = 0; j < m; ++j)
            tmp[j] = A[j * p + i];

        A_sub[i] = FFT_DivideAndConquer(tmp, ep, _1);
    }

    T ek = _1;
    for (int k = 0; k < n; ++k)
    {
        B[k] = T(0);
        T eki = _1;
        for (int i = 0; i < p; ++i)
        {
            B[k] += A_sub[i][k % m] * eki;
            eki *= ek;
        }
        ek *= e;
    }

    return B;
}

/*
    In-place FFT
*/

template <class iter, class T>
void FFT(iter first, iter last, const T e, const T _1 = T(1))
{
    const int n = std::distance(first, last);
    assert(__builtin_popcount(n) == 1);

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
    !!! n must be a power of 2 and e must be and n-root of unity (_1)
*/
template <class T>
std::vector<T> FFT_Iterative(const std::vector<T>& A,
                             const T e,
                             const T _1 = T(1))
{
    std::vector<T> B(A);
    FFT(B.begin(), B.end(), e, _1);
    return B;
}
