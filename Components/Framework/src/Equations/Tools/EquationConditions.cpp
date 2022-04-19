/**
 * @file EquationConditions.cpp
 * @brief Source of the equation condition functors
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/EquationConditions.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Conditions {

   bool IsPrognostic::operator()(std::shared_ptr<EquationData> eqA)
   {
      return eqA->equationType() == CouplingInformation::PROGNOSTIC;
   }

   bool IsDiagnostic::operator()(std::shared_ptr<EquationData> eqA)
   {
      return eqA->equationType() == CouplingInformation::DIAGNOSTIC;
   }

   bool IsTrivial::operator()(std::shared_ptr<EquationData> eqA)
   {
      return eqA->equationType() == CouplingInformation::TRIVIAL;
   }

   bool IsWrapper::operator()(std::shared_ptr<EquationData> eqA)
   {
      return eqA->equationType() == CouplingInformation::WRAPPER;
   }
}
}
}
