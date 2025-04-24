/**
 * @file OpGrouped.hpp
 * @brief Wrapper for transpose implementation backends
 */
#pragma once

// External includes
//

// Project includes
//
#include "ViewOps/Transpose/Cpu/OpGrouped.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "ViewOps/Transpose/Cuda/OpGrouped.hpp"
#endif
#ifdef QUICC_MPI
#include "ViewOps/Transpose/Mpi/OpGrouped.hpp"
#endif
