/**
 * @file P_Zero.cpp
 * @brief Source of the implementation of the Worland P_Zero integrator
 */

// System includes
//
#include <cassert>

// Debug includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Integrator/P_Zero.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   P_Zero::P_Zero()
   {
      this->mProfileId = Debug::Profiler::WORLANDINTG_P_ZERO;
   }

   P_Zero::~P_Zero()
   {
   }

   void P_Zero::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = 0; // operator shifts l by one
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

}
}
}
}
}
