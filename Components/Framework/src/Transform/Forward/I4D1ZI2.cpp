/** 
 * @file I4D1ZI2.cpp
 * @brief Source of forward projection operator I4D1ZI2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4D1ZI2.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4D1ZI2::sTag()
   {
      return "Fwd::I4D1ZI2";
   }

   const std::size_t& I4D1ZI2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4D1ZI2>(I4D1ZI2::sTag());
      return *i;
   }

   I4D1ZI2::I4D1ZI2()
      : IOperator(I4D1ZI2::sTag())
   {
   }

   I4D1ZI2::~I4D1ZI2()
   {
   }

}
}
}
