/** 
 * @file D1ZP.cpp
 * @brief Source of backward projection operator D1ZP 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1ZP.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1ZP::sTag()
   {
      return "Bwd::D1ZP";
   }

   const std::size_t& D1ZP::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D1ZP>(D1ZP::sTag());
      return *i;
   }

   D1ZP::D1ZP()
      : IOperator(D1ZP::sTag())
   {
   }

   D1ZP::~D1ZP()
   {
   }

}
}
}
