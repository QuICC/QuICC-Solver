/**
 * @file ApplyGalerkinStencil.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_APPLYGALERKINSTENCIL_HPP
#define QUICC_EQUATIONS_APPLYGALERKINSTENCIL_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/Utility.hpp"
#include "Types/Typedefs.hpp"
#include "Arithmetics/Utility.hpp"
#include "Arithmetics/LinearAlgebra.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Apply the galerkin stencil operator
    *
    * @param eq         Equation
    * @param compId     Equation field component ID
    * @param rField     Output field
    * @param start      Start index in linear storage
    * @param matIdx     System index
    * @param rhs        RHS field data
    */
   template <typename TData> void applyGalerkinStencil(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const TData& rhs);

   template <typename TData> inline void applyGalerkinStencil(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const TData& rhs)
   {
      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      auto outBlk = std::make_tuple(0, 0, op->rows(), Arithmetics::getCols(rhs));
      auto inBlk = std::make_tuple(start, 0, op->cols(), Arithmetics::getCols(rhs));
      Arithmetics::computeAx<Arithmetics::Operation::Set>(rField, outBlk, *op, rhs, inBlk);
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_APPLYGALERKINSTENCIL_HPP
