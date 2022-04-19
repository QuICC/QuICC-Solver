/** 
 * @file D1.cpp
 * @brief Source of forward projection operator D1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/D1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string D1::sTag()
   {
      return "Fwd::D1";
   }

   const std::size_t& D1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D1>(D1::sTag());
      return *i;
   }

   D1::D1()
      : IOperator(D1::sTag())
   {
   }

   D1::~D1()
   {
   }

}
}
}
