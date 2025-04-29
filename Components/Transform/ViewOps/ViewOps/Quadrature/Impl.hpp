/**
 * @file Impl.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// System includes
//

// Project includes
//
#include "ViewOps/Quadrature/Cpu/Impl.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Quadrature/Cuda/Impl.hpp"
#endif
