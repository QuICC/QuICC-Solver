/**
 * @file I4Q.cpp
 * @brief Source of forward projection operator I4Q
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4Q.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4Q::sTag()
   {
      return "Fwd::I4Q";
   }

   const std::size_t& I4Q::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4Q>(I4Q::sTag());
      return *i;
   }

   I4Q::I4Q()
      : IOperator(I4Q::sTag())
   {
   }

   I4Q::~I4Q()
   {
   }

}
}
}
