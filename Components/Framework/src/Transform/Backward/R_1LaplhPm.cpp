/** 
 * @file R_1LaplhPm.cpp
 * @brief Source of backward projection operator R_1LaplhPm 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/R_1LaplhPm.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string R_1LaplhPm::sTag()
   {
      return "Bwd::R_1LaplhPm";
   }

   const std::size_t& R_1LaplhPm::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<R_1LaplhPm>(R_1LaplhPm::sTag());
      return *i;
   }

   R_1LaplhPm::R_1LaplhPm()
      : IOperator(R_1LaplhPm::sTag())
   {
   }

   R_1LaplhPm::~R_1LaplhPm()
   {
   }

}
}
}
