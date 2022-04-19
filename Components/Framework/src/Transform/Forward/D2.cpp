/** 
 * @file D2.cpp
 * @brief Source of forward projection operator D2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/D2.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string D2::sTag()
   {
      return "Fwd::D2";
   }

   const std::size_t& D2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D2>(D2::sTag());
      return *i;
   }

   D2::D2()
      : IOperator(D2::sTag())
   {
   }

   D2::~D2()
   {
   }

}
}
}
