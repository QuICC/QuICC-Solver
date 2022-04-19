/**
 * @file D1Laplh.cpp
 * @brief Source of backward projection operator D1Laplh
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/D1Laplh.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string D1Laplh::sTag()
   {
      return "Bwd::D1Laplh";
   }

   const std::size_t& D1Laplh::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<D1Laplh>(D1Laplh::sTag());
      return *i;
   }

   D1Laplh::D1Laplh()
      : IOperator(D1Laplh::sTag())
   {
   }

   D1Laplh::~D1Laplh()
   {
   }

}
}
}
