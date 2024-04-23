/**
 * @file EnergyR2.cpp
 * @brief Source of the implementation of the Bessel R^2 energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/EnergyR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   EnergyR2<base_t>::EnergyR2()
      : EnergyReductor<PowerR2>()
   {
      this->setProfileTag();
   }

}
}
}
}
}
