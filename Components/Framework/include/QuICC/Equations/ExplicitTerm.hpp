/**
 * @file ExplicitTerm.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_EXPLICITTERM_HPP
#define QUICC_EQUATIONS_EXPLICITTERM_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Utility.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctor.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctorSSR.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctorSMR.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctorM.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctorS.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param explicitField Explicit linear field values
    * @param matIdx        System index
    * @param useShift   Use galerkin shifts
    * @param shiftTop   shift top?
    */
   template <typename T, typename TOperator, typename TData>
      void addExplicitTerm(const Resolution& res, const CouplingInformation& cinfo, const TOperator& op, TData& rSolverField, const int eqStart, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx, const bool useShift, const bool shiftTop);


   template <typename T, typename TOperator, typename TData>
      void addExplicitTerm(const Resolution& res, const CouplingInformation& cinfo, const TOperator& op, TData& rSolverField, const int eqStart, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx, const bool useShift, const bool shiftTop)
   {
      if constexpr((std::is_same_v<T,MHDFloat> || std::is_same_v<T, MHDComplex> ) && (std::is_same_v<TOperator, SparseMatrixZ> && std::is_same_v<typename Arithmetics::GetScalarType<TData>::ScalarType, MHDFloat>))
      {
         throw std::logic_error("This should not be called");
      } else
      {
         // Create pointer to sparse operator
         if(cinfo.indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
         {
            details::ExplicitTermFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(res, cinfo, matIdx, useShift, shiftTop);
            func.apply(rSolverField, op, eqStart, explicitField);
         }
         else if(cinfo.indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
         {
            details::ExplicitTermFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(res, cinfo, matIdx, useShift, shiftTop);
            func.apply(rSolverField, op, eqStart, explicitField);
         }
         else if(cinfo.indexType() == CouplingIndexType::MODE)
         {
            details::ExplicitTermFunctor<CouplingIndexType::MODE> func(res, cinfo, matIdx, useShift, shiftTop);
            func.apply(rSolverField, op, eqStart, explicitField);
         }
         else if(cinfo.indexType() == CouplingIndexType::SINGLE)
         {
            details::ExplicitTermFunctor<CouplingIndexType::SINGLE> func(res, cinfo, matIdx, useShift, shiftTop);
            func.apply(rSolverField, op, eqStart, explicitField);
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EXPLICITTERM_HPP
