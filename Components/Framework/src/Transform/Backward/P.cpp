/**
 * @file P.cpp
 * @brief Source of backward projection operator P
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/P.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string P::sTag()
   {
      return "Bwd::P";
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
