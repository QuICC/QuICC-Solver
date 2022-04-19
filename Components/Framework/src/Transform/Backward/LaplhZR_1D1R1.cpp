/** 
 * @file LaplhZR_1D1R1.cpp
 * @brief Source of backward projection operator LaplhZR_1D1R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/LaplhZR_1D1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string LaplhZR_1D1R1::sTag()
   {
      return "Bwd::LaplhZR_1D1R1";
   }

   const std::size_t& LaplhZR_1D1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<LaplhZR_1D1R1>(LaplhZR_1D1R1::sTag());
      return *i;
   }

   LaplhZR_1D1R1::LaplhZR_1D1R1()
      : IOperator(LaplhZR_1D1R1::sTag())
   {
   }

   LaplhZR_1D1R1::~LaplhZR_1D1R1()
   {
   }

}
}
}
