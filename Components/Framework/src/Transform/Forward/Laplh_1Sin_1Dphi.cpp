/** 
 * @file Laplh_1Sin_1Dphi.cpp
 * @brief Source of forward projection operator Laplh_1Sin_1Dphi 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1Sin_1Dphi.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1Sin_1Dphi::sTag()
   {
      return "Fwd::Laplh_1Sin_1Dphi";
   }

   const std::size_t& Laplh_1Sin_1Dphi::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Laplh_1Sin_1Dphi>(Laplh_1Sin_1Dphi::sTag());
      return *i;
   }

   Laplh_1Sin_1Dphi::Laplh_1Sin_1Dphi()
      : IOperator(Laplh_1Sin_1Dphi::sTag())
   {
   }

   Laplh_1Sin_1Dphi::~Laplh_1Sin_1Dphi()
   {
   }

}
}
}
