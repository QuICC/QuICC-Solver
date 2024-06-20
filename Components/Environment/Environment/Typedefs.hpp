/**
 * @file Typedefs.hpp
 * @brief Typedefs for environment
 */

#ifndef QUICC_ENVIRONMENT_TYPEDEFS_HPP
#define QUICC_ENVIRONMENT_TYPEDEFS_HPP

#ifdef QUICC_MPI
#include "mpi.h"
#endif //QUICC_MPI

namespace QuICC {

namespace Environment {

#ifdef QUICC_MPI
   typedef MPI_Comm CommType;
#else
   typedef int CommType;
#endif

} // Environment
} // QuICC
#endif // QUICC_ENVIRONMENT_TYPEDEFS_HPP
