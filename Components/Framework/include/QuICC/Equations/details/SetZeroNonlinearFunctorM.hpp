/**
 * @file SetZeroNonlinearFunctorM.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTORM_HPP
#define QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTORM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/details/SetZeroNonlinearFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void SetZeroNonlinearFunctor<CouplingIndexType::MODE>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   // Safety assertion
   assert(start >= 0);

   // Get mode indexes
   ArrayI mode = tRes.mode(matIdx);
   int rows = field.comp(compId).slice(mode(0)).rows();
   int zeroRow = cinfo.galerkinShift(matIdx, 0);

   // Set data to zero
   int k = start;
   for (int i = zeroRow; i < rows; i++)
   {
      // Set field to zero
      Arithmetics::setZero(storage, k);

      // increase storage counter
      k++;
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTORM_HPP
