/** 
 * @file Sin_1Dphi.cpp
 * @brief Source of backward projection operator Sin_1Dphi 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1Dphi.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1Dphi::sTag()
   {
      return "Bwd::Sin_1Dphi";
   }

   const std::size_t& Sin_1Dphi::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Sin_1Dphi>(Sin_1Dphi::sTag());
      return *i;
   }

   Sin_1Dphi::Sin_1Dphi()
      : IOperator(Sin_1Dphi::sTag())
   {
   }

   Sin_1Dphi::~Sin_1Dphi()
   {
   }

}
}
}
