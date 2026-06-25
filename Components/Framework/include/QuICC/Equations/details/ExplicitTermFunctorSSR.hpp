/**
 * @file ExplicitTermFunctorSSR.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORSSR_HPP
#define QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORSSR_HPP

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
void ExplicitTermFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::apply(
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
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp(op.cols(),
         1);
#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
   static_assert(false, "Parallel MPI solve is not supported anymore");
#else
      int k = 0;
      const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
      for (int j = 0; j < cols; j++)
      {
         // Effective rows in case of non-uniform truncation
         int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

         for (int i = 0; i < rows; i++)
         {
            // Copy slice into flat array
            tmp(k, 0) = explicitField.point(i, j, matIdx);

            // increase storage counter
            k++;
         }
      }
#endif // defined QUICC_MPI && defined QUICC_MPISPSOLVE

      if(zeroRow > 0 || shiftMaxRow > 0)
      {
         // Apply operator to field
         std::tuple<int, int, int, int> outBlk =
            std::make_tuple(0, 0, op.rows(), Arithmetics::getCols(tmp));
         Arithmetics::computeAx<Arithmetics::Operation::Set>(tmp, outBlk,
               op, tmp);

         int k = 0;
         int kk = eqStart;
         // effective cols
         const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
         const int usedCols = cols - shiftMaxCol;
         for (int j = 0; j < cols; j++)
         {
            // Effective rows
            int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
            int usedRows = rows - shiftMaxRow;

            for (int i = 0; i < rows; i++)
            {
               if(j >= zeroCol && i >= zeroRow && j < usedCols && i < usedRows)
               {
                  // Add data to solver field
                  Arithmetics::assignScalar<Arithmetics::Operation::Plus>(rSolverField, kk,
                        tmp(k,0));
                  kk++;
               }

               // increase storage counter
               k++;
            }
         }
      }
      else
      {
         // Apply operator to field
         std::tuple<int, int, int, int> outBlk =
            std::make_tuple(eqStart, 0, op.rows(), Arithmetics::getCols(tmp));
         Arithmetics::computeAx<Arithmetics::Operation::Plus>(rSolverField, outBlk,
               op, tmp);
      }
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORSSR_HPP
