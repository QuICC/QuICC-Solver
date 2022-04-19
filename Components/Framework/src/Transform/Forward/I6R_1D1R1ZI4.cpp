/** 
 * @file I6R_1D1R1ZI4.cpp
 * @brief Source of forward projection operator I6R_1D1R1ZI4 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6R_1D1R1ZI4.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6R_1D1R1ZI4::sTag()
   {
      return "Fwd::I6R_1D1R1ZI4";
   }

   const std::size_t& I6R_1D1R1ZI4::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I6R_1D1R1ZI4>(I6R_1D1R1ZI4::sTag());
      return *i;
   }

   I6R_1D1R1ZI4::I6R_1D1R1ZI4()
      : IOperator(I6R_1D1R1ZI4::sTag())
   {
   }

   I6R_1D1R1ZI4::~I6R_1D1R1ZI4()
   {
   }

}
}
}
