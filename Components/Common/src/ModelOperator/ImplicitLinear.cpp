/**
 * @file ImplicitLinear.cpp
 * @brief Source of the ImplicitLinear ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/ImplicitLinear.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string ImplicitLinear::sTag()
   {
      return "implicit_linear";
   }

   std::string ImplicitLinear::sFormatted()
   {
      return "ImplicitLinear";
   }

   ImplicitLinear::ImplicitLinear()
      : IRegisterId<ImplicitLinear>(ImplicitLinear::sTag(), ImplicitLinear::sFormatted())
   {
   }

}
}
