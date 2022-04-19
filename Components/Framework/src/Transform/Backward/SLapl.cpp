/**
 * @file SLapl.cpp
 * @brief Source of backward projection operator SLapl
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/SLapl.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string SLapl::sTag()
   {
      return "Bwd::SLapl";
   }

   const std::size_t& SLapl::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<SLapl>(SLapl::sTag());
      return *i;
   }

   SLapl::SLapl()
      : IOperator(SLapl::sTag())
   {
   }

   SLapl::~SLapl()
   {
   }

}
}
}
