/**
 * @file TruncateFieldFunctorSMR.hpp
 * @brief Implementation of the CopyUnknown functor fpr SLOWEST_MULTI_RHS
 */

#ifndef QUICC_EQUATIONS_DETAILS_TRUNCATEFIELDFUNCTORSMR_HPP
#define QUICC_EQUATIONS_DETAILS_TRUNCATEFIELDFUNCTORSMR_HPP

// System includes
//
#include <memory>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/details/TruncateFieldFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <bool IsSet, typename TData, typename TField>
void TruncateFieldFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);

   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx) - shiftMaxCol;

   // Safety assertion
   assert(start >= 0);

   // Copy data
   for (int j = zeroCol; j < cols; j++)
   {
      const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx) - shiftMaxRow;
      for (int i = zeroRow; i < rows; i++)
      {
         if constexpr (IsSet)
         {
            // Copy field value into storage
            Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage,
               i - zeroRow + start, j - zeroCol,
               field.point(i, j, matIdx));
         }
         else
         {
            // Add field value to storage
            Arithmetics::assignScalar<Arithmetics::Operation::Plus>(storage,
               i - zeroRow + start, j - zeroCol,
               field.point(i, j, matIdx));
         }
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_TRUNCATEFIELDFUNCTORSMR_HPP
