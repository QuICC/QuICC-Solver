/** 
 * @file P0.cpp
 * @brief Source of forward projection operator P0 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/P0.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string P0::sTag()
   {
      return "Fwd::P0";
   }

   const std::size_t& P0::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<P0>(P0::sTag());
      return *i;
   }

   P0::P0()
      : IOperator(P0::sTag())
   {
   }

   P0::~P0()
   {
   }

}
}
}
