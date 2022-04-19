/** 
 * @file I4D1.cpp
 * @brief Source of forward projection operator I4D1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4D1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4D1::sTag()
   {
      return "Fwd::I4D1";
   }

   const std::size_t& I4D1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4D1>(I4D1::sTag());
      return *i;
   }

   I4D1::I4D1()
      : IOperator(I4D1::sTag())
   {
   }

   I4D1::~I4D1()
   {
   }

}
}
}
