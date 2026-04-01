/**
 * @file StoreSolutionFunctorSMR.hpp
 * @brief Implementation of the StoreSolution functor for SLOWEST_MULTI_RHS
 */

#ifndef QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORSMR_HPP
#define QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORSMR_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void StoreSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(
   TField& field, const TData& storage, const int start)
{
   int solStart;
   auto solution = init(solStart, storage, start);

   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);

   // Copy data
   for (int j = 0; j < cols; j++)
   {
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
      for (int i = 0; i < rows; i++)
      {
         // Copy timestep output into field
         MHDVariant dataPoint =
            Arithmetics::getScalar(*solution.ptr, i + solStart, j);
         dataPoint = (*spUp)(dataPoint, i, j, matIdx);
         field.setPoint(dataPoint, i, j, matIdx);
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_STORESOLUTIONFUNCTORSMR_HPP
