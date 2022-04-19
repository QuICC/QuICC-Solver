/**
 * @file R_1.cpp
 * @brief Source of backward projection operator R_1
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1::sTag()
   {
      return "Bwd::R_1";
   }

   const std::size_t& R_1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<R_1>(R_1::sTag());
      return *i;
   }

   R_1::R_1()
      : IOperator(R_1::sTag())
   {
   }

   R_1::~R_1()
   {
   }

}
}
}
