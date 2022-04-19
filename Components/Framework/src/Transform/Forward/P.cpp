/** 
 * @file P.cpp
 * @brief Source of forward projection operator P 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/P.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string P::sTag()
   {
      return "Fwd::P";
   }

   const std::size_t& P::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<P>(P::sTag());
      return *i;
   }

   P::P()
      : IOperator(P::sTag())
   {
   }

   P::~P()
   {
   }

}
}
}
