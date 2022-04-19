/** 
 * @file D3.cpp
 * @brief Source of backward projection operator D3 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D3.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D3::sTag()
   {
      return "Bwd::D3";
   }

   const std::size_t& D3::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D3>(D3::sTag());
      return *i;
   }

   D3::D3()
      : IOperator(D3::sTag())
   {
   }

   D3::~D3()
   {
   }

}
}
}
