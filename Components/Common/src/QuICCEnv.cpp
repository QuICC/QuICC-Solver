/** 
 * @file IEnvironment.cpp
 * @brief Source of generic environment
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/QuICCEnv.hpp"

// Project includes
//
#ifdef QUICC_MPI
#include "QuICC/Environment/Mpi.hpp"
#define QUICC_ENVIMPL Environment::Mpi
#else
#include "QuICC/Environment/Serial.hpp"
#define QUICC_ENVIMPL Environment::Serial
#endif // QUICC_MPI

namespace QuICC {

   Environment::IEnvironment& QuICCEnv()
   {
      static QUICC_ENVIMPL env;
      return env;
   }

#undef QUICC_ENVIMPL

}
