/** 
 * @file D2.cpp
 * @brief Source of backward projection operator D2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D2.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D2::sTag()
   {
      return "Bwd::D2";
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
