/**
 * @file EquationSorters.cpp
 * @brief Source of the equation sorting functors
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/EquationSorters.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Sorters {

   bool EquationType::operator()(std::shared_ptr<EquationData> eqA, std::shared_ptr<EquationData> eqB)
   {
      return computeEquationType(eqA) < computeEquationType(eqB);
   }

   int EquationType::computeEquationType(std::shared_ptr<EquationData> eqA)
   {
      return static_cast<int>(eqA->equationType());
   }
}
}
}
