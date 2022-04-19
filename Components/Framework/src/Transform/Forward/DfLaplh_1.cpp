/** 
 * @file DfLaplh_1.cpp
 * @brief Source of forward projection operator DfLaplh_1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/DfLaplh_1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string DfLaplh_1::sTag()
   {
      return "Fwd::DfLaplh_1";
   }

   const std::size_t& DfLaplh_1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<DfLaplh_1>(DfLaplh_1::sTag());
      return *i;
   }

   DfLaplh_1::DfLaplh_1()
      : IOperator(DfLaplh_1::sTag())
   {
   }

   DfLaplh_1::~DfLaplh_1()
   {
   }

}
}
}
