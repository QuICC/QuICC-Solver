/** 
 * @file R_2.cpp
 * @brief Source of backward projection operator R_2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_2.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_2::sTag()
   {
      return "Bwd::R_2";
   }

   const std::size_t& R_2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<R_2>(R_2::sTag());
      return *i;
   }

   R_2::R_2()
      : IOperator(R_2::sTag())
   {
   }

   R_2::~R_2()
   {
   }

}
}
}
