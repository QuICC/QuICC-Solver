/**
 * @file Stencil.cpp
 * @brief Source of the Stencil ModelOperator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperator/Stencil.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperator {

   std::string Stencil::sTag()
   {
      return "stencil";
   }

   std::string Stencil::sFormatted()
   {
      return "Stencil";
   }

   Stencil::Stencil()
      : IRegisterId<Stencil>(Stencil::sTag(), Stencil::sFormatted())
   {
   }

}
}
