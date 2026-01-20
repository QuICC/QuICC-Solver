/**
 * @file AddSourceFunctorSMR.hpp
 * @brief AddSourceFunctor implementation For SLOWEST_MULTI_RHS
 */

#ifndef QUICC_EQUATIONS_ADDSOURCEFUNCTORSMR_HPP
#define QUICC_EQUATIONS_ADDSOURCEFUNCTORSMR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/AddSourceFunctor.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

   template <>
   template <typename TData, typename TField> void AddSourceFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& info = eq->couplingInfo(compId);
      // Add source term if required
      if(info.hasSource())
      {
         const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
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

         //Safety assertion
         assert(start >= 0);

         // Copy data
         for(int j = zeroCol; j < cols; j++)
         {
            int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
            for(int i = zeroRow; i < rows; i++)
            {
               // Add source term
               Arithmetics::assignScalar<Arithmetics::Operation::Plus>(storage, i - zeroRow + start, j - zeroCol, eq->sourceTerm(compId, i, j, matIdx));
            }
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCEFUNCTORSMR_HPP
