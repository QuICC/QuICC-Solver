/** 
 * @file Sin_1D1Sin.cpp
 * @brief Source of backward projection operator Sin_1D1Sin 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1D1Sin.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1D1Sin::sTag()
   {
      return "Bwd::Sin_1D1Sin";
   }

   const std::size_t& Sin_1D1Sin::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Sin_1D1Sin>(Sin_1D1Sin::sTag());
      return *i;
   }

   Sin_1D1Sin::Sin_1D1Sin()
      : IOperator(Sin_1D1Sin::sTag())
   {
   }

   Sin_1D1Sin::~Sin_1D1Sin()
   {
   }

}
}
}
