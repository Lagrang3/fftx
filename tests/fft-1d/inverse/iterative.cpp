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

typedef complex<double> cd;
const double PI = acos(-1.0);

int main()
{
    int n, N = 1;
    cin >> n;

    N = n;
    // while (N < n) N <<= 1;

    vector<cd> A(N, 0);

    for (int i = 0; i < n; ++i)
    {
        double x;
        cin >> x;
        A[i] = cd(x, 0);
    }
    auto FA = FFT_Iterative(A, cd(cos(2 * PI / N), sin(2 * PI / N)));
    auto A2 = FFT_Iterative(FA, cd(cos(2 * PI / N), -sin(2 * PI / N)));
    double inv_n = 1.0 / N;
    for (auto& x : A2)
        x *= inv_n;
    // done

    // get differences
    double diff = 0;
    for (int i = 0; i < N; ++i)
    {
        diff += norm(A2[i] - A[i]);
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
