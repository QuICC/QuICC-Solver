/**
 * @file SetBoundaryValueFunctorSMR.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTORSMR_HPP
#define QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTORSMR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/SetBoundaryValueFunctor.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

template <>
   template <typename TData, typename TField> void SetBoundaryValueFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      const auto& info = eq->couplingInfo(compId);
      const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
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
         const int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
         for(int i = zeroRow; i < rows; i++)
         {
            // Add source term
            Arithmetics::assignScalar<Arithmetics::Operation::Set>(storage, i - zeroRow + start, j - zeroCol, eq->boundaryValue(compId, i, j, matIdx));
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTORSMR_HPP
