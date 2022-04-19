/** 
 * @file I6Laplh.cpp
 * @brief Source of forward projection operator I6Laplh 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/I6Laplh.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string I6Laplh::sTag()
   {
      return "Fwd::I6Laplh";
   }

   const std::size_t& I6Laplh::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<I6Laplh>(I6Laplh::sTag());
      return *i;
   }

   I6Laplh::I6Laplh()
      : IOperator(I6Laplh::sTag())
   {
   }

   I6Laplh::~I6Laplh()
   {
   }

}
}
}
