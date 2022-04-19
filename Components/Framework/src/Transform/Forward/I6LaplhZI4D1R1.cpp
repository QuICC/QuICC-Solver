/** 
 * @file I6LaplhZI4D1R1.cpp
 * @brief Source of forward projection operator I6LaplhZI4D1R1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6LaplhZI4D1R1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6LaplhZI4D1R1::sTag()
   {
      return "Fwd::I6LaplhZI4D1R1";
   }

   const std::size_t& I6LaplhZI4D1R1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I6LaplhZI4D1R1>(I6LaplhZI4D1R1::sTag());
      return *i;
   }

   I6LaplhZI4D1R1::I6LaplhZI4D1R1()
      : IOperator(I6LaplhZI4D1R1::sTag())
   {
   }

   I6LaplhZI4D1R1::~I6LaplhZI4D1R1()
   {
   }

}
}
}
