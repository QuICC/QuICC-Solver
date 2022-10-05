/**
 * @file EnergyD1R1.cpp
 * @brief Source of the reductor transform operator Reductor::EnergyD1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergyD1R1::sTag()
   {
      return "Red::EnergyD1R1";
   }

   std::string EnergyD1R1::sFormatted()
   {
      return "Reductor::EnergyD1R1";
   }

   EnergyD1R1::EnergyD1R1()
      : IRegisterId<EnergyD1R1>(EnergyD1R1::sTag(), EnergyD1R1::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
