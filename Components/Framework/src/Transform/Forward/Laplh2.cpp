/** 
 * @file Laplh2.cpp
 * @brief Source of forward projection operator Laplh2 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh2.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh2::sTag()
   {
      return "Fwd::Laplh2";
   }

   const std::size_t& Laplh2::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Laplh2>(Laplh2::sTag());
      return *i;
   }

   Laplh2::Laplh2()
      : IOperator(Laplh2::sTag())
   {
   }

   Laplh2::~Laplh2()
   {
   }

}
}
}
