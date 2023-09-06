/**
 * @file R1_Zero.cpp
 * @brief Source of the implementation of the Worland R1_Zero integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/Base/R1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   R1_Zero<base_t>::R1_Zero()
   {
      this->setProfileTag();
   }


   void R1_Zero<base_t>::initBackend() const
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
