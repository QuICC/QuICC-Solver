/**
 * @file ExplicitTermFunctorM.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_EXPLICITTERMFUNCTORM_HPP
#define QUICC_EQUATIONS_EXPLICITTERMFUNCTORM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/ExplicitTermFunctor.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

template <>
   template <typename T, typename TOperator, typename TData>
   void ExplicitTermFunctor<CouplingIndexType::MODE>::compute(const std::size_t opId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField)
   {
      if constexpr((std::is_same<T,MHDFloat>::value || std::is_same<T, MHDComplex>::value ) && (std::is_same<TOperator, SparseMatrixZ>::value && std::is_same<TData, Matrix>::value))
      {
      } else
      {
         // Create pointer to sparse operator
         const TOperator * op = &eq->template explicitOperator<TOperator>(opId, compId, fieldId, matIdx);

         const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         // Get mode indexes
         ArrayI mode = tRes.mode(matIdx);

         // Assert correct sizes
         assert(op->cols() == explicitField.slice(mode(0)).rows());

         // Apply operator to field
         Arithmetics::addMatrixProduct(rSolverField, eqStart, *op, explicitField.slice(mode(0)).col(mode(1)));
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EXPLICITTERMFUNCTORM_HPP
