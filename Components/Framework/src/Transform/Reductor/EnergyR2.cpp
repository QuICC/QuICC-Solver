/** 
 * @file EnergyR2.cpp
 * @brief Source of reductor projection operator EnergyR2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergyR2.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergyR2::sTag()
   {
      return "Red::EnergyR2";
   }

   const std::size_t& EnergyR2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<EnergyR2>(EnergyR2::sTag());
      return *i;
   }

   EnergyR2::EnergyR2()
      : IOperator(EnergyR2::sTag())
   {
   }

   EnergyR2::~EnergyR2()
   {
   }

}
}
}
