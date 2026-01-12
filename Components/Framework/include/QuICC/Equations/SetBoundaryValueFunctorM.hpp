/**
 * @file SetBoundaryValueFunctorM.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTORM_HPP
#define QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTORM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/SetBoundaryValueFunctor.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

template <>
   template <typename TData, typename TField> void SetBoundaryValueFunctor<CouplingIndexType::MODE>::apply(const TField& field, TData& storage, const int start)
   {
      const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      //Safety assertion
      assert(start >= 0);

      // Get mode indexes
      ArrayI mode = tRes.mode(matIdx);
      int rows = field.comp(compId).slice(mode(0)).rows();
      int zeroRow = eq->couplingInfo(compId).galerkinShift(matIdx,0);

      // Copy data
      int k = start;
      for(int i = zeroRow; i < rows; i++)
      {
         // Add source term
         Arithmetics::setScalar(storage, k, eq->boundaryValue(compId, i, mode(1), mode(0)));

         // increase storage counter
         k++;
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTORM_HPP
