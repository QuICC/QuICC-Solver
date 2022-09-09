/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Worland energy operator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Worland/Reductor/Energy.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   Energy::Energy()
      : EnergyReductor<Power>()
   {
      this->setProfileTag();
   }

   Energy::~Energy()
   {
   }

}
}
}
}
}
