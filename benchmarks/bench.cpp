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

double test_fftw(const vector<cd>& data, size_t N, int R)
{
    double T = 0;

    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    for (size_t i = 0; i < N; ++i)
    {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int r = 0; r < R; ++r)
    {
        auto t1 = high_resolution_clock::now();
        fftw_execute(p);
        auto dt = high_resolution_clock::now() - t1;
        T += duration_cast<microseconds>(dt).count();
    }
    fftw_free(in);
    fftw_free(out);
    fftw_destroy_plan(p);

    return T / R / 1000;
}

template <class Ftype>
double test_mylib(const vector<cd>& data, size_t N, int R, Ftype F)
{
    double T = 0;
    const cd e(cos(2 * PI / N), -sin(2 * PI / N));
    vector<cd> A(N);
    copy(data.begin(), data.begin() + N, A.begin());

    for (int r = 0; r < R; ++r)
    {
        auto t1 = high_resolution_clock::now();
        F(A, e, cd(1));
        auto dt = high_resolution_clock::now() - t1;
        T += duration_cast<microseconds>(dt).count();
    }
    /*
    template <class T>
    std::vector<T> FFT_BruteForce(const std::vector<T>& A,
                                  const T e,
                                  const T _1 = T(1))
    template <class T>
    std::vector<T> FFT_DivideAndConquer(const std::vector<T>& A,
                                        const T e,
                                        const T _1 = T(1))
    template <class T>
    std::vector<T> FFT_Iterative(const std::vector<T>& A,
                                 const T e,
                                 const T _1 = T(1))
    */
    return T / R / 1000;
}

int main()
{
    const int Nmax = 1 << 20;
    uniform_real_distribution<double> U(0, 1);
    vector<cd> data(Nmax);

    for (auto& x : data)
        x = cd(U(rng), U(rng));

    for (int N = 1 << 8, R; N < 2000; N *= 4)
    {
        R = 100;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "FFTW: " << test_fftw(data, N, R)
             << " ms\n";
        cout << setw(30) << left << "MyFFT BruteForce: "
             << test_mylib(data, N, R, FFT_BruteForce<cd>) << " ms\n";
        cout << setw(30) << left << "MyFFT DivideAndConquer: "
             << test_mylib(data, N, R, FFT_DivideAndConquer<cd>) << " ms\n";
        cout << setw(30) << left
             << "MyFFT InPlace: " << test_mylib(data, N, R, FFT_InPlace<cd>)
             << " ms\n";
        cout << setw(30) << left
             << "MyFFT Iterative: " << test_mylib(data, N, R, FFT_Iterative<cd>)
             << " ms\n";
        cout << "\n\n";
    }
    for (int N = 1 << 14, R; N <= Nmax; N *= 4)
    {
        R = 10;
        cout << "N = " << N << '\n';
        cout << setw(30) << left << "FFTW: " << test_fftw(data, N, R)
             << " ms\n";
        cout << setw(30) << left << "MyFFT DivideAndConquer: "
             << test_mylib(data, N, R, FFT_DivideAndConquer<cd>) << " ms\n";
        cout << setw(30) << left
             << "MyFFT InPlace: " << test_mylib(data, N, R, FFT_InPlace<cd>)
             << " ms\n";
        cout << setw(30) << left
             << "MyFFT Iterative: " << test_mylib(data, N, R, FFT_Iterative<cd>)
             << " ms\n";
        cout << "\n\n";
    }

    return 0;
}
