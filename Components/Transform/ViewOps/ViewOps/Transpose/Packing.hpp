/**
 * @file Packing.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// External includes
//

// Project includes
//
// #include "ViewOps/Transpose/Cpu/Packing.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Transpose/Cuda/Packing.hpp"
#endif
