/**
 * @file Op.hpp
 * @brief Wrapper for implementation backends
 */
#pragma once

// System includes
//

// Project includes
//
#include "ViewOps/Slicewise/Cpu/NoGridOp.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
//#include "ViewOps/Slicewise/Cuda/NoGridOp.hpp"
#endif
