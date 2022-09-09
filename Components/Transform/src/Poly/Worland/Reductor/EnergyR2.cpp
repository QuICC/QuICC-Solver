/**
 * @file EnergyR2.cpp
 * @brief Source of the implementation of the Worland R^2 energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   EnergyR2::EnergyR2()
      : EnergyReductor<PowerR2>()
   {
      this->setProfileTag();
   }

   EnergyR2::~EnergyR2()
   {
   }

}
}
}
}
}
