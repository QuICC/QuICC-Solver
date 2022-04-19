/** 
 * @file Mpi.cpp
 * @brief Source of MPI initializer and finalizer
 */

// Configuration includes
//

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "QuICC/Environment/Register/Mpi.hpp"

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Environment/IEnvironment.hpp"

namespace QuICC {

namespace Environment {

namespace Register {

   const int Mpi::priority = 10;

   const int Mpi::initId = IEnvironment::registerInitializer(Mpi::priority,&Mpi::initializer);

   const int Mpi::finalId = IEnvironment::registerFinalizer(Mpi::priority,&Mpi::finalizer);

   void Mpi::initializer()
   {
      DebuggerMacro_msg("Running registered MPI library initializer", 1);

      MPI_Init(0, 0);
   }

   void Mpi::finalizer()
   {
      DebuggerMacro_msg("Running registered MPI library finalizer", 1);

      MPI_Finalize();
   }
}
}
}
