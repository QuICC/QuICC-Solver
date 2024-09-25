/**
 * @file Op.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// System includes
//

// Project includes
//
#include "ViewOps/Slicewise/Cpu/Op.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Slicewise/Cuda/Op.hpp"
#endif
