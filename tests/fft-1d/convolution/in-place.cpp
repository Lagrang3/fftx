/*
    Test fft library with complex<double>,
    the test consist in the finding the convolution
    of two arrays A and B using brute force and
    using FFT; then we compare differences.

    Eduardo Quintana Miranda
    31-05-2019
*/

#include <cmath>
#include <complex>
#include <fftx.hpp>
#include <iostream>
using namespace std;
using namespace fftx;

typedef complex<double> cd;
const double PI = acos(-1.0);

int main()
{
    int n, N = 1;
    cin >> n;

    while (N < n)
        N <<= 1;
    N <<= 1;

    vector<cd> A(N, 0), B(N, 0), C(N);

    for (int i = 0; i < n; ++i)
    {
        double x;
        cin >> x;
        A[i] = cd(x, 0);
    }
    for (int i = 0; i < n; ++i)
    {
        double x;
        cin >> x;
        B[i] = cd(x, 0);
    }

    // brute force convolution
    for (int i = 0; i < N; ++i)
    {
        C[i] = 0;
        for (int j = 0; j <= i; ++j)
            C[i] += A[j] * B[i - j];
    }
    // done

    // FFT-based convolution
    auto FA = FFT_InPlace(A, cd(cos(2 * PI / N), sin(2 * PI / N)));
    auto FB = FFT_InPlace(B, cd(cos(2 * PI / N), sin(2 * PI / N)));

    vector<cd> FC(N, 0);
    for (int i = 0; i < N; ++i)
        FC[i] = FA[i] * FB[i];

    auto C2 = FFT_InPlace(FC, cd(cos(2 * PI / N), -sin(2 * PI / N)));
    double inv_n = 1.0 / N;
    for (auto& x : C2)
        x *= inv_n;
    // done

    // get differences
    double diff = 0;
    for (int i = 0; i < N; ++i)
    {
        diff += norm(C2[i] - C[i]);
    }
    diff = sqrt(diff) / N;
    // done

    const double eps = 1e-8;

    if (diff > eps)
    {
        cout << "Final coefficients differ by: " << diff << "\n";
        return 1;
    }

    return 0;
}
