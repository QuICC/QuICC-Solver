/**
 * @file Reduction.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// System includes
//

// Project includes
//
#include "ViewOps/Reduction/Cpu/Reduction.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Reduction/Cuda/Reduction.hpp"
#endif
