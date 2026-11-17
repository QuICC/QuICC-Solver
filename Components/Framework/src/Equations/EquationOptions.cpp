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
      : EquationOptions(0, true, true, true)
   {
   }

   EquationOptions::EquationOptions(const int it)
      : EquationOptions(it, true, true, true)
   {
   }

   EquationOptions::EquationOptions(const int it, const bool nlIsLhs, const bool traHasQi)
      : EquationOptions(it, nlIsLhs, traHasQi, true)
   {
   }

   EquationOptions::EquationOptions(const int it, const bool nlIsLhs, const bool traHasQi, const bool isBase)
      : nonlinearIsLhs(nlIsLhs), transformHasQi(traHasQi), isBase(isBase), mIt(it)
   {
   }

   int EquationOptions::it() const
   {
      return this->mIt;
   }

} // Equations
} // QuICC
