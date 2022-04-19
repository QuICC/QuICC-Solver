/** 
 * @file Hdf5.cpp
 * @brief Source of HDF5 initializer and finalizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Environment/Register/Hdf5.hpp"

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Environment/IEnvironment.hpp"

namespace QuICC {

namespace Environment {

namespace Register {

   const int Hdf5::priority = 20;

   const int Hdf5::initId = IEnvironment::registerInitializer(Hdf5::priority,&Hdf5::initializer);

   const int Hdf5::finalId = IEnvironment::registerFinalizer(Hdf5::priority,&Hdf5::finalizer);

   void Hdf5::initializer()
   {
      DebuggerMacro_msg("Running registered HDF5 library initializer", 1);
   }

   void Hdf5::finalizer()
   {
      DebuggerMacro_msg("Running registered HDF5 library finalizer", 1);

      H5close();
   }
}
}
}
