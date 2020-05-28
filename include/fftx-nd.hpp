#pragma once

#include <fftx-1d.hpp>
#include <memory>

namespace fftx
{
template <int dim, class iter, class T>
void FFT_dim(iter first, iter last, const int N, const T e, const T _1 = T(1))
{
    size_t len = power(N, dim);
    std::vector<T> buff(len);

    std::copy(first, last, buff.begin());

    for (int d = 0; d < dim; ++d)
    {
        // perform 1D fft
        for (size_t i = 0; i < buff.size();)
        {
            int j = i + N;
            FFT_InPlace(&buff[i], &buff[j], e, _1);
            i = j;
        }

        const int pN = power(N, dim - 1);
        // transpose
        std::vector<T> tmp(len);
        for (size_t i = 0, j; i < tmp.size(); ++i)
        {
            j = i / N + (i % N) * pN;
            tmp[i] = buff[j];
        }

        buff = std::move(tmp);
    }
    std::copy(buff.begin(), buff.end(), first);
}
}
