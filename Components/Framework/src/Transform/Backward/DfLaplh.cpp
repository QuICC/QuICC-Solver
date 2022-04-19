/** 
 * @file DfLaplh.cpp
 * @brief Source of backward projection operator DfLaplh 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/DfLaplh.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string DfLaplh::sTag()
   {
      return "Bwd::DfLaplh";
   }

   const std::size_t& DfLaplh::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<DfLaplh>(DfLaplh::sTag());
      return *i;
   }

   DfLaplh::DfLaplh()
      : IOperator(DfLaplh::sTag())
   {
   }

   DfLaplh::~DfLaplh()
   {
   }

}
}
}
