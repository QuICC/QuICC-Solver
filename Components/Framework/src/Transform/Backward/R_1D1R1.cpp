/**
 * @file R_1D1R1.cpp
 * @brief Source of backward projection operator R_1D1R1
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1D1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1D1R1::sTag()
   {
      return "Bwd::R_1D1R1";
   }

   const std::size_t& R_1D1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<R_1D1R1>(R_1D1R1::sTag());
      return *i;
   }

   R_1D1R1::R_1D1R1()
      : IOperator(R_1D1R1::sTag())
   {
   }

   R_1D1R1::~R_1D1R1()
   {
   }

}
}
}
