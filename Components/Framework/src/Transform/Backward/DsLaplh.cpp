/** 
 * @file DsLaplh.cpp
 * @brief Source of backward projection operator DsLaplh 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/DsLaplh.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string DsLaplh::sTag()
   {
      return "Bwd::DsLaplh";
   }

   const std::size_t& DsLaplh::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<DsLaplh>(DsLaplh::sTag());
      return *i;
   }

   DsLaplh::DsLaplh()
      : IOperator(DsLaplh::sTag())
   {
   }

   DsLaplh::~DsLaplh()
   {
   }

}
}
}
