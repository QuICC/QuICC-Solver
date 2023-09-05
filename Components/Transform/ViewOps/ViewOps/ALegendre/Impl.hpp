/**
 * @file Impl.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// External includes
//

// Project includes
//
#include "ViewOps/ALegendre/Cpu/Impl.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/ALegendre/Cuda/Impl.hpp"
#endif
// #include "ViewOps/ALegendre/Tags.hpp"
