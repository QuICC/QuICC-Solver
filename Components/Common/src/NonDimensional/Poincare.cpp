/**
 * @file Poincare.cpp
 * @brief Source of the Poincare nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Poincare.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Poincare::sTag()
   {
      return "poincare";
   }

   std::string Poincare::sFormatted()
   {
      return "Poincare";
   }

   Poincare::Poincare(const MHDFloat value)
      : IRegisterId<Poincare>(value, Poincare::sTag(), Poincare::sFormatted())
   {
   }

}
}
