/**
 * @file QuICCEnv.cpp
 * @brief Source of static environment
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Environment/QuICCEnv.hpp"

// Project includes
//
#ifdef QUICC_MPI
#include "Environment/Mpi.hpp"
#define QUICC_ENVIMPL Environment::Mpi
#else
#include "Environment/Serial.hpp"
#define QUICC_ENVIMPL Environment::Serial
#endif // QUICC_MPI

namespace QuICC {

   /**
    * @brief Enviroment singleton
    */
   Environment::IEnvironment& QuICCEnv()
   {
      static QUICC_ENVIMPL env;

      return env;
   }

#undef QUICC_ENVIMPL

}
