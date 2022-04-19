/** 
 * @file Sin_1Laplh.cpp
 * @brief Source of backward projection operator Sin_1Laplh 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Sin_1Laplh.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Sin_1Laplh::sTag()
   {
      return "Bwd::Sin_1Laplh";
   }

   const std::size_t& Sin_1Laplh::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Sin_1Laplh>(Sin_1Laplh::sTag());
      return *i;
   }

   Sin_1Laplh::Sin_1Laplh()
      : IOperator(Sin_1Laplh::sTag())
   {
   }

   Sin_1Laplh::~Sin_1Laplh()
   {
   }

}
}
}
