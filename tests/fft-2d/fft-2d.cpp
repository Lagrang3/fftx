#include <cmath>
#include <complex>
#include <fftx.hpp>
#include <iostream>
#include <random>

#include <fftw3.h>

using namespace std;

typedef complex<double> cd;
const double PI = acos(-1.0);
default_random_engine rng;

int main()
{
    for (int N = 1; N < 2000; N <<= 1)
    {
        vector<cd> A(N * N, 0);
        const cd e(cos(2 * PI / N), -sin(2 * PI / N));
        uniform_real_distribution<double> U(0, 1);

        for (auto& x : A)
            x = cd(U(rng), U(rng));

        // my fft-2d
        vector<cd> FA(A);
        FFT_dim<2>(FA.begin(), FA.end(), N, e);

        // fftw3-2d
        fftw_complex *in, *out;
        fftw_plan p;
        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * A.size());
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * A.size());

        for (size_t i = 0; i < A.size(); ++i)
        {
            in[i][0] = A[i].real();
            in[i][1] = A[i].imag();
        }
        p = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);

        double diff = 0, eps = 1e-10 * FA.size();
        for (size_t i = 0; i < FA.size(); ++i)
        {
            diff += abs(FA[i] - cd(out[i][0], out[i][1]));
        }

        fftw_free(in);
        fftw_free(out);
        fftw_destroy_plan(p);

        if (diff > eps)
            return 1;
    }
    return 0;
}
