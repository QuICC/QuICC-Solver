/**
 * @file Op.hpp
 * @brief Wrapper for transpose implementation backends
 */
#pragma once

// External includes
//

// Project includes
//
#include "ViewOps/Transpose/Cpu/Op.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Transpose/Cuda/Op.hpp"
#endif
#ifdef QUICC_MPI
#include "ViewOps/Transpose/Mpi/Op.hpp"
#endif
