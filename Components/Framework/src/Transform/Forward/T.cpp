/** 
 * @file T.cpp
 * @brief Source of forward projection operator T 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/T.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string T::sTag()
   {
      return "Fwd::T";
   }

   const std::size_t& T::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<T>(T::sTag());
      return *i;
   }

   T::T()
      : IOperator(T::sTag())
   {
   }

   T::~T()
   {
   }

}
}
}
