/**
 * @file Transpose.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// External includes
//

// Project includes
//
#include "ViewOps/Transpose/Cpu/Transpose.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Transpose/Cuda/Transpose.hpp"
#endif
#ifdef QUICC_MPI
#include "ViewOps/Transpose/Mpi/Transpose.hpp"
#endif

