/**
 * @file ViewMacros.hpp
 * @brief Common macros for View library
 */

#pragma once

// System includes
//

// Project includes
//

#ifdef __CUDACC__
#define QUICC_CUDA_HOSTDEV __host__ __device__
#define QUICC_CUDA_HOST __host__
#else
#define QUICC_CUDA_HOSTDEV
#define QUICC_CUDA_HOST
#endif
