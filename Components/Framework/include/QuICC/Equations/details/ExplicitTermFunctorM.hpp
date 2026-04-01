/**
 * @file ExplicitTermFunctorM.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORM_HPP
#define QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORM_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/LinearAlgebra.hpp"
#include "Arithmetics/Utility.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctor.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename T, typename TOperator, typename TData>
void ExplicitTermFunctor<CouplingIndexType::MODE>::apply(
   TData& rSolverField, const TOperator& op, const int eqStart,
   const typename Framework::Selector::ScalarField<T>& explicitField)
{
   if constexpr ((std::is_same<T, MHDFloat>::value ||
                    std::is_same<T, MHDComplex>::value) &&
                 (std::is_same<TOperator, SparseMatrixZ>::value &&
                    std::is_same<TData, Matrix>::value))
   {}
   else
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      // Get mode indexes
      ArrayI mode = tRes.mode(matIdx);

      // Assert correct sizes
      assert(op.cols() == explicitField.slice(mode(0)).rows());

      // Apply operator to field
      typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp =
         explicitField.slice(mode(0)).col(mode(1)).eval();
      std::tuple<int, int, int, int> outBlk =
         std::make_tuple(eqStart, 0, op.rows(), Arithmetics::getCols(tmp));
      Arithmetics::computeAx<Arithmetics::Operation::Plus>(rSolverField, outBlk,
         op, tmp);
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORM_HPP
