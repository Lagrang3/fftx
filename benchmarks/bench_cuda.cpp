#include <chrono>
#include <cmath>
#include <complex>
#include <fftx.hpp>
#include <iomanip>
#include <iostream>
#include <random>

#include <fftw3.h>

using namespace std;
using namespace std::chrono;

typedef complex<double> cd;
const double PI = acos(-1.0);
default_random_engine rng;
const int Nmax = 1 << 24;
uniform_real_distribution<double> U(0, 1);
vector<cd> _data(Nmax);


std::vector<cd> cuFFT(const std::vector<cd>& data,const int N);

double test_cuda(const vector<cd>& data, size_t N,const int R)
{
    double T = 0;

    for (int r = 0; r < R; ++r)
    {
        auto t1 = high_resolution_clock::now();
        cuFFT(data,N);
        auto dt = high_resolution_clock::now() - t1;
        T += duration_cast<microseconds>(dt).count();
    }
    return T / R / 1000;
}

void benchmark_radix2()
{
    cout << "    << Radix-2 >>"
         << "\n";
    for (int N = 1 << 8, R; N < 2000; N *= 4)
    {
        R = 2;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "cuFFT: " << test_cuda(_data, N, R)
             << " ms\n";
        cout << "\n\n";
    }
    for (int N = 1 << 14, R; N <= Nmax; N *= 4)
    {
        R = 2;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "cuFFT: " << test_cuda(_data, N, R)
             << " ms\n";
        cout << "\n\n";
    }
}
void benchmark_any_radix()
{
    cout << "    << any Radix >>"
         << "\n";
    for (int N = 100, R; N < 2000; N *= 10)
    {
        R = 2;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "cuFFT: " << test_cuda(_data, N,R)
             << " ms\n";
        cout << "\n\n";
    }
    for (int N = 10000, R; N <= Nmax; N *= 10)
    {
        R = 2;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "FFTW: " << test_cuda(_data, N,R)
             << " ms\n";
        cout << "\n\n";
    }
}
void benchmark_primes()
{
    vector<int> primes{109, 211, 401, 809, 1009, 10009};
    cout << "    << Primes >>"
         << "\n";
    for (auto N : primes)
    {
        int R = 2;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "cuFFT: " << test_cuda(_data, N, R)
             << " ms\n";
        cout << "\n\n";
    }
}

int main()
{
    for (auto& x : _data)
        x = cd(U(rng), U(rng));

    benchmark_radix2();
    benchmark_any_radix();
    benchmark_primes();

    return 0;
}
