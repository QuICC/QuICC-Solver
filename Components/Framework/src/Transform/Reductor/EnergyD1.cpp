/** 
 * @file EnergyD1.cpp
 * @brief Source of reductor projection operator EnergyD1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergyD1.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergyD1::sTag()
   {
      return "Red::EnergyD1";
   }

   const std::size_t& EnergyD1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<EnergyD1>(EnergyD1::sTag());
      return *i;
   }

   EnergyD1::EnergyD1()
      : IOperator(EnergyD1::sTag())
   {
   }

   EnergyD1::~EnergyD1()
   {
   }

}
}
}
