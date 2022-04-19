/** 
 * @file SRadLapl.cpp
 * @brief Source of backward projection operator SRadLapl 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/SRadLapl.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string SRadLapl::sTag()
   {
      return "Bwd::SRadLapl";
   }

   const std::size_t& SRadLapl::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<SRadLapl>(SRadLapl::sTag());
      return *i;
   }

   SRadLapl::SRadLapl()
      : IOperator(SRadLapl::sTag())
   {
   }

   SRadLapl::~SRadLapl()
   {
   }

}
}
}
