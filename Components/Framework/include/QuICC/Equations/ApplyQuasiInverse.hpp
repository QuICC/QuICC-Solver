/**
 * @file ApplyQuasiInverse.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_APPLYQUASIINVERSE_HPP
#define QUICC_EQUATIONS_APPLYQUASIINVERSE_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Arithmetics/Basic.hpp"
#include "Arithmetics/LinearAlgebra.hpp"
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Apply the quasi-inverse operator
    *
    * @param rField     Output field
    * @param start      Start index in linear storage
    * @param matIdx     System index
    * @param rhsStart   Start index in RHS data
    * @param rhs        RHS field data
    */
   template <typename TOperator, typename TData> void applyQuasiInverse(const TOperator& op, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs, const bool isSet = false);
   template <bool IsSet, typename TOp, typename TData> void applyQuasiInverseImpl(const TOp& op, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs);

   template <typename TOperator, typename TData> inline void applyQuasiInverse(const TOperator& op, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs, const bool isSet)
   {
      if(isSet)
      {
         applyQuasiInverseImpl<true>(op, rField, start, matIdx, rhsStart, rhs);
      }
      else
      {
         applyQuasiInverseImpl<false>(op, rField, start, matIdx, rhsStart, rhs);
      }
   }

   template <bool IsSet, typename TOp, typename TData> inline void applyQuasiInverseImpl(const TOp& op, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs)
   {
      // Get number of rows and cols
      int cols = Arithmetics::getCols(rField);
      int rhsRows = op.cols();

      auto outBlk = std::make_tuple(start, 0, op.rows(), Arithmetics::getCols(rhs));
      auto inBlk = std::make_tuple(rhsStart, 0, rhsRows, cols);
      if constexpr(IsSet)
      {
         Arithmetics::computeAx<Arithmetics::Operation::Set>(rField, outBlk, op, rhs, inBlk);
      }
      else
      {
         Arithmetics::computeAx<Arithmetics::Operation::Plus>(rField, outBlk, op, rhs, inBlk);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_APPLYQUASIINVERSE_HPP
