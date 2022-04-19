/** 
 * @file GpuDctTools.hpp
 * @brief Utility GPU kernels for DCT transform
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUDCTTOOLS_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUDCTTOOLS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

namespace Gpu {

   /**
    * @brief Build input for DCT as DFT
    */
   __global__ void buildDCTInput(cufftDoubleReal* out, const cufftDoubleReal* in, const int size, const int blockSize);

   /**
    * @brief Compute DCT phases
    */
   __global__ void buildDCTPhases(double* s, double* cs, const int fftSize, const int size);

   /**
    * @brief Compute DCT phases
    */
   __global__ void buildDCT4Phases(double* s, double* cs, const int fftSize, const int size);

   /**
    * @brief Extract DCT output from DFT
    */
   __global__ void buildDCTOutput(cufftDoubleReal* out, const cufftDoubleComplex* in, const int size, const int blockSize);

   /**
    * @brief Extract DCT output from DFT
    */
   __global__ void buildDCTOutput(cufftDoubleReal* out, const cufftDoubleComplex* in, const int size, const int blockSize, const double * s, const double* cs);

   /**
    * @brief Build input for IDCT as DFT
    */
   __global__ void buildIDCTInput(cufftDoubleComplex* out, const cufftDoubleReal* in, const int size, const int blockSize);

   /**
    * @brief Build input for IDCT as DFT
    */
   __global__ void buildIDCTInput(cufftDoubleComplex* out, const cufftDoubleReal* in, const int size, const int blockSize, const double *sn, const double *cs);

   /**
    * @brief Extract IDCT output from DFT
    */
   __global__ void buildIDCTOutput(cufftDoubleReal* out, const cufftDoubleReal* in, const int size, const int blockSize);

   /**
    * @brief Build input for DCT-IV as DFT
    */
   __global__ void buildDCT4Input(cufftDoubleComplex* out, const cufftDoubleReal* in, const int size, const int blockSize);

   /**
    * @brief Build input for DCT-IV as DFT
    */
   __global__ void buildDCT4Input(cufftDoubleComplex* out, const cufftDoubleReal* in, const int size, const int blockSize, const double *sn, const double *cs);

   /**
    * @brief Extract DCT-IV output from DFT
    */
   __global__ void buildDCT4Output(cufftDoubleReal* out, const cufftDoubleComplex* in, const int size, const int blockSize, const double scale);

   /**
    * @brief Extract DCT-IV output from DFT
    */
   __global__ void buildDCT4Output(cufftDoubleReal* out, const cufftDoubleComplex* in, const int size, const int blockSize, const double scale, const double *sn, const double *cs);

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_GPUDCTTOOLS_HPP
