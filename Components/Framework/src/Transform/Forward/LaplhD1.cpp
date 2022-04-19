/** 
 * @file LaplhD1.cpp
 * @brief Source of forward projection operator LaplhD1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/LaplhD1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string LaplhD1::sTag()
   {
      return "Fwd::LaplhD1";
   }

   const std::size_t& LaplhD1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<LaplhD1>(LaplhD1::sTag());
      return *i;
   }

   LaplhD1::LaplhD1()
      : IOperator(LaplhD1::sTag())
   {
   }

   LaplhD1::~LaplhD1()
   {
   }

}
}
}
