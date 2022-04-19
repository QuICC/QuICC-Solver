/** 
 * @file Laplh_1.cpp
 * @brief Source of forward projection operator Laplh_1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1::sTag()
   {
      return "Fwd::Laplh_1";
   }

   const std::size_t& Laplh_1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Laplh_1>(Laplh_1::sTag());
      return *i;
   }

   Laplh_1::Laplh_1()
      : IOperator(Laplh_1::sTag())
   {
   }

   Laplh_1::~Laplh_1()
   {
   }

}
}
}
