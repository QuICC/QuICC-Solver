/**
 * @file I4S.cpp
 * @brief Source of forward projection operator I4S
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4S.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4S::sTag()
   {
      return "Fwd::I4S";
   }

   const std::size_t& I4S::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4S>(I4S::sTag());
      return *i;
   }

   I4S::I4S()
      : IOperator(I4S::sTag())
   {
   }

   I4S::~I4S()
   {
   }

}
}
}
