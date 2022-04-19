/** 
 * @file Sin_1.cpp
 * @brief Source of backward projection operator Sin_1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1::sTag()
   {
      return "Bwd::Sin_1";
   }

   const std::size_t& Sin_1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Sin_1>(Sin_1::sTag());
      return *i;
   }

   Sin_1::Sin_1()
      : IOperator(Sin_1::sTag())
   {
   }

   Sin_1::~Sin_1()
   {
   }

}
}
}
