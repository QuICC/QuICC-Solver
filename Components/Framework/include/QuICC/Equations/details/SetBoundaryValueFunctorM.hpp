/**
 * @file SetBoundaryValueFunctorM.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORM_HPP
#define QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORM_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void SetBoundaryValueFunctor<CouplingIndexType::MODE>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   // Safety assertion
   assert(start >= 0);

   // Get mode indexes
   ArrayI mode = tRes.mode(matIdx);
   int rows = field.slice(mode(0)).rows();
   int zeroRow = cinfo.galerkinShift(matIdx, 0);

   // Copy data
   int k = start;
   for (int i = zeroRow; i < rows; i++)
   {
      // Add source term
      Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage, k,
         spBoundary->compute(i, mode(1), mode(0)));

      // increase storage counter
      k++;
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTORM_HPP
