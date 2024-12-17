/**
 * @file CudaUtil.hpp
 * @brief
 */

#pragma once

// External includes
//
#include <stdexcept>
#include <cuda_runtime_api.h>


// Project includes
//

namespace QuICC {
/// @brief This namespace provides cuda specific utilities
namespace Cuda {

/// Macro to enable error checking
#ifdef NDEBUG
#   define cudaErrChk(ans) { ans; }
#else
#   define cudaErrChk(ans) { QuICC::Cuda::cudaAssert((ans), __FILE__, __LINE__); }
#endif
/// @brief Utility to check errors
/// @param cErr cuda error returned by the invoked function
/// @param file location of invoking function
/// @param line location of invoking function
inline void cudaAssert(cudaError_t cErr, const char *file, int line)
{
    if (cErr)
    {
        constexpr unsigned int size = 200;
        char msg [size];
        snprintf(msg, size,
          "cudaAssert: %s: %s:%d\n", cudaGetErrorString(cErr), file, line);
        throw  std::runtime_error(msg);
    }
}

/// @brief Check if the memory is allocated on the device
/// @param ptr pointer to memory location
/// @return true if device, false if host
inline bool isDeviceMemory(const void* ptr)
{
  cudaPointerAttributes att;
  cudaErrChk(cudaPointerGetAttributes(&att, ptr));
  return att.type == cudaMemoryTypeDevice;
}


} // namespace Cuda
} // namespace QuICC

