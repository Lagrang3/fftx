#include <stdexcept>
#include <thrust/device_vector.h>
#include <vector>
#include <complex>
#include <cufft.h>

using cd = std::complex<double>;

auto cuFFT(const std::vector<cd>& data,const int N)
{
    
    
    thrust::device_vector<cd> datain;
    
    thrust::copy(data.begin(),data.begin()+N, std::back_inserter(datain) );
    
    cufftDoubleComplex * dataptr = reinterpret_cast<cufftDoubleComplex*>(
        thrust::raw_pointer_cast(&datain[0]));
    cufftHandle plan;
    
    if(cufftPlan1d(&plan,datain.size(),CUFFT_Z2Z,1)!=CUFFT_SUCCESS)
    {
        throw std::runtime_error("Cuda error: plan creation failed");
    }
    
    if(cufftExecZ2Z(plan,dataptr,dataptr,CUFFT_FORWARD)!=CUFFT_SUCCESS)
    {
        throw std::runtime_error("Cuda error: ExecZ2Z failed");
    }
    
    if(cudaThreadSynchronize()!=cudaSuccess)
    {   
        throw std::runtime_error("Cuda error: failed to syncronize");
    }
    cufftDestroy(plan);
    
    std::vector<cd> dataout;
    thrust::copy(datain.begin(),datain.end(), std::back_inserter(dataout) );
    return dataout;
}
