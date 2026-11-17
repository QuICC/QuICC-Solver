/**
 * @file CorrectSolutionFunctorM.cpp
 * @brief Implementation of the CorrectSolution functor fpr MODE
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
void CorrectSolutionFunctor<CouplingIndexType::MODE>::init(const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections)
{
   for (auto&& c: corrections)
   {
      const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      ArrayI mode = tRes.mode(matIdx);
      if (mode(1) == std::get<2>(c) && mode(0) == std::get<3>(c))
      {
         corr.emplace_back(std::get<MHDComplex>(std::get<0>(c)), std::get<1>(c), 0);
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC
