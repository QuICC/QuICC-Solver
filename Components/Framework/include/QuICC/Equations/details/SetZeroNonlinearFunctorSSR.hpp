/**
 * @file SetZeroNonlinearFunctorSSR.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTORSSR_HPP
#define QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTORSSR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Arithmetics/Basic.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/details/SetZeroNonlinearFunctor.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename TData, typename TField>
void SetZeroNonlinearFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(
   const TField& field, TData& storage, const int start)
{
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   // Safety assertion
   assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   for (int k = 0; k < info.galerkinN(matIdx); ++k)
   {
      // Set field to zero
      Arithmetics::setZero(storage, k + start);
   }
#else
   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
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

   // Set data to zero
   int k = start;
   for (int j = zeroCol; j < cols; j++)
   {
      // Effective rows in case of non-uniform truncation
      int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

      for (int i = zeroRow; i < usedRows; i++)
      {
         // Set field to zero
         Arithmetics::setZero(storage, k);

         // increase storage counter
         k++;
      }
   }
#endif // defined QUICC_MPI && defined QUICC_MPISPSOLVE
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTORSSR_HPP
