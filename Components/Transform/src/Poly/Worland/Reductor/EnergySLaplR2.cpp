/**
 * @file EnergySLaplR2.cpp
 * @brief Source of the implementation of the Worland Spherical Laplacian R^2 energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/EnergySLaplR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   EnergySLaplR2::EnergySLaplR2()
      : EnergyReductor<PowerSLaplR2>()
   {
      this->setProfileTag();
   }

   EnergySLaplR2::~EnergySLaplR2()
   {
   }

}
}
}
}
}
