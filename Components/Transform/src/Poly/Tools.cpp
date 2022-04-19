/**
 * @file Tools.cpp
 * @brief Defines some useful constants and tools for polynomial transforms
 */

// System includes
//
#include <cmath>

// External includes
//

// Class include
//
#include "QuICC/Transform/Poly/Tools.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Poly {

   const MHDFloat Tools::STD_DEALIASING = 3.0/2.0;

   int Tools::dealias(const int size)
   {
      return std::ceil(Tools::STD_DEALIASING*static_cast<MHDFloat>(size));
   }

}
}
}
