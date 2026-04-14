/**
 * @file EquationOptions.cpp
 * @brief Source of base class to hold equation options
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/EquationOptions.hpp"

namespace QuICC {

namespace Equations {

   EquationOptions::EquationOptions()
      : EquationOptions(0, true)
   {
   }

   EquationOptions::EquationOptions(const int it)
      : EquationOptions(it, true)
   {
   }

   EquationOptions::EquationOptions(const int it, const bool nlIsLhs)
      : nonlinearIsLhs(nlIsLhs), mIt(it)
   {
   }

   int EquationOptions::it() const
   {
      return this->mIt;
   }

} // Equations
} // QuICC
