/**
 * @file EnergyD1.cpp
 * @brief Source of the reductor transform operator Reductor::EnergyD1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergyD1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergyD1::sTag()
   {
      return "Red::EnergyD1";
   }

   std::string EnergyD1::sFormatted()
   {
      return "Reductor::EnergyD1";
   }

   EnergyD1::EnergyD1()
      : IRegisterId<EnergyD1>(EnergyD1::sTag(), EnergyD1::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
