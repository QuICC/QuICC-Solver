/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Worland energy operator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Energy.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

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
