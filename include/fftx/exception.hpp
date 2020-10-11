#pragma once

#include <exception>

namespace fftx
{
    class error : public std::runtime_error
    {
       public:
        using std::runtime_error ::runtime_error;
    };
}  // namespace fftx
