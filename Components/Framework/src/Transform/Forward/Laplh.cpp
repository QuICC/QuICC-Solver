/** 
 * @file Laplh.cpp
 * @brief Source of forward projection operator Laplh 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh::sTag()
   {
      return "Fwd::Laplh";
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
