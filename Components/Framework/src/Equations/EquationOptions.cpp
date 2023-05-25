/**
 * @file EquationOptions.cpp
 * @brief Source of base class to hold equation options
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Equations/EquationOptions.hpp"

namespace QuICC {

namespace Equations {

   EquationOptions::EquationOptions(const int it)
      : mIt(it)
   {
   }

   int EquationOptions::it() const
   {
      return this->mIt;
   }

} // Equations
} // QuICC
