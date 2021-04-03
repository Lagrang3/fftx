#pragma once

/*
    Modular Arithmetics
*/

#include <iostream>

namespace my_modulo_lib
{
    template <class T, class Exp>
    T power(const T& x, Exp n, const T& I = T{1})
    {
        T r{I};
        for (T aux = x; n; n /= 2)
        {
            if (n % 2)
                r *= aux;
            aux *= aux;
        }
        return r;
    }
    template <typename T, T x>
    class field_modulo
    {
       public:
        typedef T integer;
        static constexpr T mod = x;
    };

    /*
        Modular Integers
    */
    template <typename Field>
    class mint;

    template <typename Field>
    std::ostream& operator<<(std::ostream& os, const mint<Field>& A)
    {
        return os << A.x << " (mod " << Field::mod << ")";
    }
    template <typename Field>
    class mint
    {
        typedef typename Field::integer integer;
        integer x;

       public:
        constexpr mint() : x{0} {}

        template <typename int_type>
        mint(int_type _x) : x{_x}
        {
            x %= Field::mod;
            if (x < 0)
                x += Field::mod;
        }

        mint(const mint& that) : x{that.x} {}

        mint& operator=(const mint& that) { return x = that.x, *this; }

        // explicit operator bool() const { return x == integer{0}; }
        operator integer() const { return x; }

        mint& operator+=(const mint& that)
        {
            return x = (x + that.x) % Field::mod, *this;
        }

        mint& operator*=(const mint& t)
        {
            // direct multiplication
            x = (x * t.x) % Field::mod;
            return *this;
        }
        bool operator==(const mint& that) const { return x == that.x; }
        bool operator!=(const mint& that) const { return x != that.x; }

        auto inverse() const { return power(*this, Field::mod - 2); }

        friend std::ostream& operator<<<Field>(std::ostream& os,
                                               const mint<Field>& A);
    };

    template <typename T>
    mint<T> operator+(const mint<T>& A, const mint<T>& B)
    {
        mint<T> C{A};
        return C += B;
    }
    template <typename T>
    mint<T> operator*(const mint<T>& A, const mint<T>& B)
    {
        mint<T> C{A};
        return C *= B;
    }
    template <typename T>
    mint<T>& operator/=(mint<T>& A, const mint<T>& B)
    {
        return A *= B.inverse();
    }
    template <typename T>
    mint<T> operator/(const mint<T>& A, const mint<T>& B)
    {
        mint<T> C{A};
        return C /= B;
    }
}  // namespace my_modulo_lib
