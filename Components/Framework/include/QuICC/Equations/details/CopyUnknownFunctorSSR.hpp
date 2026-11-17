/**
 * @file CopyUnknownFunctorSSR.hpp
 * @brief Implementation of the CopyUnknown functor for SLOWEST_SINGLE_RHS
 */

#ifndef QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTORSSR_HPP
#define QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTORSSR_HPP

// System includes
//
#include <memory>
#include <stdexcept>
#include <vector>

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/SetZeroNonlinear.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <bool IsSet, typename TData, typename TField>
void CopyUnknownFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);

   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx) - shiftMaxCol;

   // Safety assertion
   assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   static_assert(false, "Parallel MPI solve is not supported anymore");
#else
   // Copy data
   int k = start;
   for (int j = zeroCol; j < cols; j++)
   {
      // Effective rows in case of non-uniform truncation
      int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx) - shiftMaxRow;

      for (int i = zeroRow; i < usedRows; i++)
      {
         if constexpr (IsSet)
         {
            // Copy field value into storage
            Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage, k,
               field.comp(compId).point(i, j, matIdx));
         }
         else
         {
            // Add field value to storage
            Arithmetics::assignScalar<Arithmetics::Operation::Plus>(storage, k,
               field.comp(compId).point(i, j, matIdx));
         }

         // increase storage counter
         k++;
      }
   }
#endif // defined QUICC_MPI && defined QUICC_MPISPSOLVE
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTORSSR_HPP
