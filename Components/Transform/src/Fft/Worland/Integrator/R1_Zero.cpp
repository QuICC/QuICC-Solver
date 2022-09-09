/**
 * @file R1_Zero.cpp
 * @brief Source of the implementation of the Worland R1_Zero integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Integrator/R1_Zero.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   R1_Zero::R1_Zero()
   {
      this->setProfileTag();
   }

   R1_Zero::~R1_Zero()
   {
   }

   void R1_Zero::initBackend() const
   {
      // l = 0 mode is set to sero
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = -1; // operator shifts l by -1
      int extraN = 1; // 1 extra modes is required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

}
}
}
}
}
