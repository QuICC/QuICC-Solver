/**
 * @file EnergyD1R1.cpp
 * @brief Source of the implementation of the Bessel D R energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/EnergyD1R1.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

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
