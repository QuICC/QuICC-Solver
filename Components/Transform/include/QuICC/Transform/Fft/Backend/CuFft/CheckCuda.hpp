/**
 * @file CheckCuda.h
 * @brief Some CUDA and CUDA libraries helpers to check status
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHECKCUDA_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHECKCUDA_HPP

// External includes
//
#include "cublas_v2.h"
#include "cusparse.h"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

void CheckCuda(cudaError_t status, const int line, const bool useException = true);

void CheckCuda(cublasStatus_t status, const int line, const bool useException = true);

void CheckCuda(cusparseStatus_t status, const int line, const bool useException = true);

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHECKCUDA_HPP
