/** 
 * @file Profiler.cpp
 * @brief Source of profiler initializer and finalizer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Environment/Register/Profiler.hpp"

// Project includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Environment/IEnvironment.hpp"

namespace QuICC {

namespace Environment {

namespace Register {

   const int Profiler::priority = 1;

   const int Profiler::initId = IEnvironment::registerInitializer(Profiler::priority,&Profiler::initializer);

   const int Profiler::finalId = IEnvironment::registerFinalizer(Profiler::priority,&Profiler::finalizer);

   void Profiler::initializer()
   {
      DebuggerMacro_msg("Running registered profiler initializer", 1);

      ProfilerMacro_init();
   }

   void Profiler::finalizer()
   {
      DebuggerMacro_msg("Running registered profiler finalizer", 1);
   }
}
}
}
