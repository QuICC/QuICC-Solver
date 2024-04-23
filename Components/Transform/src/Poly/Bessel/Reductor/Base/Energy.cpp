/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Bessel energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Energy.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   Energy<base_t>::Energy()
      : EnergyReductor<Power>()
   {
      this->setProfileTag();
   }

}
}
}
}
}
