/** 
 * @file I4R_1D1R1.cpp
 * @brief Source of forward projection operator I4R_1D1R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I4R_1D1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I4R_1D1R1::sTag()
   {
      return "Fwd::I4R_1D1R1";
   }

   const std::size_t& I4R_1D1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I4R_1D1R1>(I4R_1D1R1::sTag());
      return *i;
   }

   I4R_1D1R1::I4R_1D1R1()
      : IOperator(I4R_1D1R1::sTag())
   {
   }

   I4R_1D1R1::~I4R_1D1R1()
   {
   }

}
}
}
