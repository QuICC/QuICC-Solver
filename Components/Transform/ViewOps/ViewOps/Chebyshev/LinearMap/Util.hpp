/**
 * @file Util.hpp
 * @brief Chebyshev LinearMap operators utilities
 */
#pragma once

// System includes
//
#include <cstdint>
#include <type_traits>

// Project includes
//
#include "ViewOps/Chebyshev/LinearMap/Tags.hpp"


#ifdef __CUDACC__
#define QUICC_CUDA_HOSTDEV __host__ __device__
#else
#define QUICC_CUDA_HOSTDEV
#endif

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {

namespace dealias {} // namespace dealias

} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
