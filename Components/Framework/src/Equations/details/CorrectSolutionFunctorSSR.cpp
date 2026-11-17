/**
 * @file CorrectSolutionFunctorSSR.cpp
 * @brief Implementation of the CorrectSolution functor fpr SLOWEST_SINGLE_RHS
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
void CorrectSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::init(const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections)
{
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);

   for (auto&& c: corrections)
   {
      if (matIdx == std::get<3>(c))
      {
         // Compute storage index
         int k = 0;
         for(int j = 0; j < std::get<2>(c); j++)
         {
            // Effective rows in case of non-uniform truncation
            k += tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
         }

         corr.emplace_back(std::get<MHDComplex>(std::get<0>(c)), k + std::get<1>(c), 0);
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC
