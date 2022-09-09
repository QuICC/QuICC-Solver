/**
 * @file EnergySLAPLR2.cpp
 * @brief Source of the implementation of the Worland Spherical Laplacian R^2 energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/EnergySLAPLR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   EnergySLAPLR2::EnergySLAPLR2()
      : EnergyReductor<PowerSLAPLR2>()
   {
      this->setProfileTag();
   }

   EnergySLAPLR2::~EnergySLAPLR2()
   {
   }

}
}
}
}
}
