/**
 * @file EnergySLaplR2.cpp
 * @brief Source of the implementation of the Worland Spherical Laplacian R^2 energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergySLaplR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   EnergySLaplR2<base_t>::EnergySLaplR2()
      : EnergyReductor<PowerSLaplR2>()
   {
      this->setProfileTag();
   }

}
}
}
}
}
