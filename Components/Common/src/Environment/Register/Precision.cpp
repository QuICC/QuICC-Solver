/** 
 * @file Precision.cpp
 * @brief Source of multiple precision initializer and finalizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Environment/Register/Precision.hpp"

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Environment/IEnvironment.hpp"
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Environment {

namespace Register {

   const int Precision::priority = 2;

   const int Precision::initId = IEnvironment::registerInitializer(Precision::priority,&Precision::initializer);

   const int Precision::finalId = IEnvironment::registerFinalizer(Precision::priority,&Precision::finalizer);

   void Precision::initializer()
   {
      DebuggerMacro_msg("Running registered multiple precision initializer", 1);

      QuICC::Precision::init();
   }

   void Precision::finalizer()
   {
      DebuggerMacro_msg("Running registered multiple precision finalizer", 1);
   }
}
}
}
