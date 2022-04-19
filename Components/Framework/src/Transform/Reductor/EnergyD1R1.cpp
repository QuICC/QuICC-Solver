/** 
 * @file EnergyD1R1.cpp
 * @brief Source of reductor projection operator EnergyD1R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Reductor/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string EnergyD1R1::sTag()
   {
      return "Red::EnergyD1R1";
   }

   const std::size_t& EnergyD1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<EnergyD1R1>(EnergyD1R1::sTag());
      return *i;
   }

   EnergyD1R1::EnergyD1R1()
      : IOperator(EnergyD1R1::sTag())
   {
   }

   EnergyD1R1::~EnergyD1R1()
   {
   }

}
}
}
