/** 
 * @file LaplhSin_1.cpp
 * @brief Source of forward projection operator LaplhSin_1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/LaplhSin_1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string LaplhSin_1::sTag()
   {
      return "Fwd::LaplhSin_1";
   }

   const std::size_t& LaplhSin_1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<LaplhSin_1>(LaplhSin_1::sTag());
      return *i;
   }

   LaplhSin_1::LaplhSin_1()
      : IOperator(LaplhSin_1::sTag())
   {
   }

   LaplhSin_1::~LaplhSin_1()
   {
   }

}
}
}
