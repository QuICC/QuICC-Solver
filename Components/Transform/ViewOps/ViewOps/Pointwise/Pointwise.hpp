/**
 * @file Pointwise.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// System includes
//

// Project includes
//
#include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Pointwise/Cuda/Pointwise.hpp"
#endif
