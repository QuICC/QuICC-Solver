/** 
 * @file LaplhSin_1Dphi.cpp
 * @brief Source of forward projection operator LaplhSin_1Dphi 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/LaplhSin_1Dphi.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string LaplhSin_1Dphi::sTag()
   {
      return "Fwd::LaplhSin_1Dphi";
   }

   const std::size_t& LaplhSin_1Dphi::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<LaplhSin_1Dphi>(LaplhSin_1Dphi::sTag());
      return *i;
   }

   LaplhSin_1Dphi::LaplhSin_1Dphi()
      : IOperator(LaplhSin_1Dphi::sTag())
   {
   }

   LaplhSin_1Dphi::~LaplhSin_1Dphi()
   {
   }

}
}
}
