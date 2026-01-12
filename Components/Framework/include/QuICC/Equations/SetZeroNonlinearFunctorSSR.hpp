/**
 * @file SetZeroNonlinearFunctorSSR.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_SETZERONONLINEARFUNCTORSSR_HPP
#define QUICC_EQUATIONS_SETZERONONLINEARFUNCTORSSR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/SetZeroNonlinearFunctor.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

template <>
template <typename TData, typename TField> void SetZeroNonlinearFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(const TField& field, TData& storage, const int start)
{
   const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
   const auto& info = eq->couplingInfo(compId);
   //Safety assertion
   assert(start >= 0);

#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   for(int k = 0; k < info.galerkinN(matIdx); ++k)
   {
      // Set field to zero
      Arithmetics::setZero(storage, k + start);
   }
#else
   int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
   int zeroRow = info.galerkinShift(matIdx,0);
   int zeroCol;
   if(eq->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      zeroCol = info.galerkinShift(matIdx,2);
   } else
   {
      zeroCol = info.galerkinShift(matIdx,1);
   }

   // Set data to zero
   int k = start;
   for(int j = zeroCol; j < cols; j++)
   {
      // Effective rows in case of non-uniform truncation
      int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

      for(int i = zeroRow; i < usedRows; i++)
      {
         // Set field to zero
         Arithmetics::setZero(storage, k);

         // increase storage counter
         k++;
      }
   }
#endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETZERONONLINEARFUNCTORSSR_HPP
