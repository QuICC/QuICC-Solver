/**
 * @file Utility.hpp
 * @brief Patch missing features from std utility
 */

#pragma once

// External includes
//
#include <type_traits>
#ifdef __CUDACC__
#include <cuda_runtime_api.h>
#endif

// Project includes
//

#ifdef __CUDACC__
#define QUICC_CUDA_HOSTDEV __host__ __device__
#define QUICC_CUDA_HOST __host__
#else
#define QUICC_CUDA_HOSTDEV
#define QUICC_CUDA_HOST
#endif


namespace cuda {
namespace std {

template<class T>
QUICC_CUDA_HOSTDEV constexpr ::std::add_const_t<T>& as_const(T& t) noexcept
{
    return t;
}


} // namespace std
} // namespace cuda

