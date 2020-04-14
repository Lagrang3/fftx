#pragma once

/*
    Fast power
    O(Log n)
*/

template <class T>
T power(const T& x, int n, const T& I = T(1))
{
    T r = I, aux = x;
    for (; n; n >>= 1)
    {
        if (n & 1)
            r *= aux;
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
