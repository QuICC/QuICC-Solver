/** 
 * @file D1ZP0.cpp
 * @brief Source of forward projection operator D1ZP0 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/D1ZP0.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string D1ZP0::sTag()
   {
      return "Fwd::D1ZP0";
   }

   const std::size_t& D1ZP0::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D1ZP0>(D1ZP0::sTag());
      return *i;
   }

   D1ZP0::D1ZP0()
      : IOperator(D1ZP0::sTag())
   {
   }

   D1ZP0::~D1ZP0()
   {
   }

}
}
}
