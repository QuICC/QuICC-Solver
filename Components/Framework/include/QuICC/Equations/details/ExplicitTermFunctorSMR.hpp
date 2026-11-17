/**
 * @file ExplicitTermFunctorSMR.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORSMR_HPP
#define QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORSMR_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Basic.hpp"
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
void ExplicitTermFunctor<CouplingIndexType::SLOWEST_MULTI_RHS>::apply(
   TData& rSolverField, const TOperator& op, const int eqStart, 
   const typename Framework::Selector::ScalarField<T>& explicitField)
{
   if constexpr ((std::is_same<T, MHDFloat>::value ||
                    std::is_same<T, MHDComplex>::value) &&
                 (std::is_same<TOperator, SparseMatrixZ>::value &&
                    std::is_same<TData, Matrix>::value))
   {}
   else if constexpr ((std::is_same<T, MHDFloat>::value) &&
                 (std::is_same<TOperator, SparseMatrixZ>::value))
   {
      throw std::logic_error("This should never be called");
   }
   else
   {
      if(zeroRow > 0 || shiftMaxRow > 0)
      {
         const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
         typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp(op.cols(),
            Arithmetics::getCols(explicitField.slice(matIdx)));

         // Apply operator to field
         std::tuple<int, int, int, int> outBlk =
            std::make_tuple(0, 0, op.rows(), Arithmetics::getCols(tmp));
         Arithmetics::computeAx<Arithmetics::Operation::Set>(tmp, outBlk,
               op, explicitField.slice(matIdx).eval());

         const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx) - shiftMaxCol;
         for (int j = zeroCol; j < cols; j++)
         {
            // Effective rows in case of non-uniform truncation
            int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx) - shiftMaxRow;

            for (int i = zeroRow; i < usedRows; i++)
            {
               // Add data to solver field
               Arithmetics::assignScalar<Arithmetics::Operation::Plus>(rSolverField, eqStart + i - zeroRow, j - zeroCol,
                     tmp(i,j));
            }
         }
      }
      else
      {
         // Apply operator to field
         std::tuple<int, int, int, int> outBlk = std::make_tuple(eqStart, 0,
            op.rows(), Arithmetics::getCols(explicitField.slice(matIdx)));
         Arithmetics::computeAx<Arithmetics::Operation::Plus>(rSolverField, outBlk,
            op, explicitField.slice(matIdx).eval());
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORSMR_HPP
