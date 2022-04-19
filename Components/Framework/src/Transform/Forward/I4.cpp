/** 
 * @file I4.cpp
 * @brief Source of forward projection operator I4 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4::sTag()
   {
      return "Fwd::I4";
   }

   const std::size_t& I4::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4>(I4::sTag());
      return *i;
   }

   I4::I4()
      : IOperator(I4::sTag())
   {
   }

   I4::~I4()
   {
   }

}
}
}
