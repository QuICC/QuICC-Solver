/**
 * @file D1R1.cpp
 * @brief Source of backward projection operator D1R1
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1R1::sTag()
   {
      return "Bwd::D1R1";
   }

   const std::size_t& D1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D1R1>(D1R1::sTag());
      return *i;
   }

   D1R1::D1R1()
      : IOperator(D1R1::sTag())
   {
   }

   D1R1::~D1R1()
   {
   }

}
}
}
