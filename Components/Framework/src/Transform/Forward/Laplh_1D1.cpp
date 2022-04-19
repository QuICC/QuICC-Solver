/** 
 * @file Laplh_1D1.cpp
 * @brief Source of forward projection operator Laplh_1D1 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/Laplh_1D1.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   std::string Laplh_1D1::sTag()
   {
      return "Fwd::Laplh_1D1";
   }

   const std::size_t& Laplh_1D1::id()
   {
      static std::size_t *i = new std::size_t();
      *i = registerId<Laplh_1D1>(Laplh_1D1::sTag());
      return *i;
   }

   Laplh_1D1::Laplh_1D1()
      : IOperator(Laplh_1D1::sTag())
   {
   }

   Laplh_1D1::~Laplh_1D1()
   {
   }

}
}
}
