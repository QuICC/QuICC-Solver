/**
 * @file EnergyR2.cpp
 * @brief Source of the implementation of the Worland R^2 energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/EnergyR2.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

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
