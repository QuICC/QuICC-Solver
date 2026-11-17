/**
 * @file CorrectSolutionFunctorSMR.cpp
 * @brief Implementation of the CorrectSolution functor fpr SLOWEST_MULTI_RHS
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/details/CorrectSolutionFunctor.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
void CorrectSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::init(const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections)
{
   for (auto&& c: corrections)
   {
      if (matIdx == std::get<3>(c))
      {
         corr.emplace_back(std::get<MHDComplex>(std::get<0>(c)), std::get<1>(c), std::get<2>(c));
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC
