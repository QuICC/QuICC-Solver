/**
 * @file P_Zero.cpp
 * @brief Source of the implementation of the Worland P_Zero integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Integrator/P_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   P_Zero<base_t>::P_Zero()
   {
      this->setProfileTag();
   }


   void P_Zero<base_t>::initBackend() const
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
