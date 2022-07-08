/** 
 * @file I4P.cpp
 * @brief Source of forward projection operator I4P 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4P.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4P::sTag()
   {
      return "Fwd::I4P";
   }

   const std::size_t& I4P::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4P>(I4P::sTag());
      return *i;
   }

   I4P::I4P()
      : IOperator(I4P::sTag())
   {
   }

   I4P::~I4P()
   {
   }

}
}
}
