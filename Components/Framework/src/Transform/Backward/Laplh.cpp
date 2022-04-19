/** 
 * @file Laplh.cpp
 * @brief Source of backward projection operator Laplh 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/Laplh.hpp"

// Project includes
//
#include "QuICC/Transform/Backward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   std::string Laplh::sTag()
   {
      return "Bwd::Laplh";
   }

   const std::size_t& Laplh::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Laplh>(Laplh::sTag());
      return *i;
   }

   Laplh::Laplh()
      : IOperator(Laplh::sTag())
   {
   }

   Laplh::~Laplh()
   {
   }

}
}
}
