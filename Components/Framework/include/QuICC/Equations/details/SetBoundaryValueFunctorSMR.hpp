/**
 * @file SetBoundaryValueFunctorSMR.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORSMR_HPP
#define QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORSMR_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void SetBoundaryValueFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
   int zeroRow = cinfo.galerkinShift(matIdx, 0);
   int zeroCol;
   if (res.sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      zeroCol = cinfo.galerkinShift(matIdx, 2);
   }
   else
   {
      zeroCol = cinfo.galerkinShift(matIdx, 1);
   }

   // Safety assertion
   assert(start >= 0);

   // Copy data
   for (int j = zeroCol; j < cols; j++)
   {
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
      for (int i = zeroRow; i < rows; i++)
      {
         // Add source term
         Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage,
            i - zeroRow + start, j - zeroCol,
            spBoundary->compute(i, j, matIdx));
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORSMR_HPP
