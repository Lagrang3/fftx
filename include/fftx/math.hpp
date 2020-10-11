#pragma once

#include <vector>

#include <fftx/exception.hpp>
/*
    Fast power
    O(Log n)
*/

namespace fftx
{
    template <class T>
    T power(const T& x, int n)
    {
        /*
            for our use case n should always be >0.
            However, if n==0 we would expect something
            like 1.
        */
        if (n < 1)
            throw fftx::error("power(x,n) expects n>1.");

        bool identity = true;
        T r{}, aux{x};

        for (; n; n >>= 1)
        {
            if (n & 1)
            {
                r = identity ? aux : r * aux;
                identity = false;
            }
            aux *= aux;
        }
        return r;
    }

    /*
        Least prime factor of n
        O(sqrt(n))
    */
    int prime_factor(int n)
    {
        if (n <= 3)
            return n;
        if (n % 2 == 0)
            return 2;
        for (int i = 3; i * i <= n; i += 2)
        {
            if (n % i == 0)
                return i;
        }
        return n;
    }

    /* Prime factorization */
    auto prime_factorization(int n)
    {
        std::vector<int> F;

        for (int x = 2; x * x <= n;)
            if (n % x == 0)
            {
                F.push_back(x);
                n /= x;
            }
            else
            {
                ++x;
            }
        if (n > 1)
            F.push_back(n);
        return F;
    }

}  // namespace fftx
