/**
 * @file EnergyD1R1.cpp
 * @brief Source of the implementation of the Worland D R energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyD1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   EnergyD1R1::EnergyD1R1()
      : EnergyReductor<PowerD1R1>()
   {
      this->setProfileTag();
   }

   EnergyD1R1::~EnergyD1R1()
   {
   }

}
}
}
}
}
