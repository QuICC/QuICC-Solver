/**
 * @file EnergyR2.cpp
 * @brief Source of the reductor transform operator Reductor::EnergyR2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergyR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergyR2::sTag()
   {
      return "Red::EnergyR2";
   }

   std::string EnergyR2::sFormatted()
   {
      return "Reductor::EnergyR2";
   }

   EnergyR2::EnergyR2()
      : IRegisterId<EnergyR2>(EnergyR2::sTag(), EnergyR2::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
