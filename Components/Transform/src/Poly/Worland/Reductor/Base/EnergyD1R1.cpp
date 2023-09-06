/**
 * @file EnergyD1R1.cpp
 * @brief Source of the implementation of the Worland D R energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   EnergyD1R1<base_t>::EnergyD1R1()
      : EnergyReductor<PowerD1R1>()
   {
      this->setProfileTag();
   }

}
}
}
}
}
