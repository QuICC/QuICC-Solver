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
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/ExplicitTermFunctor.hpp"
#include "QuICC/Equations/ExplicitTermFunctorSSR.hpp"
#include "QuICC/Equations/ExplicitTermFunctorSMR.hpp"
#include "QuICC/Equations/ExplicitTermFunctorM.hpp"
#include "QuICC/Equations/ExplicitTermFunctorS.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param eq            Equation
    * @param compId        Equation field component ID
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param fieldId       Physical field ID
    * @param explicitField Explicit linear field values
    * @param matIdx        System index
    */
   template <typename T, typename TData>
      void addExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx);
   template <typename T, typename TOperator,typename TData> void computeExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Framework::Selector::ScalarField<T>& explicitField, const int matIdx);


   template <typename T, typename TData>
      void addExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZTerm(opId, compId, fieldId))
      {
         computeExplicitTerm<T,SparseMatrixZ>(eq, opId, compId, rSolverField,  eqStart, fieldId, explicitField, matIdx);
      }

      // Compute with real linear operator
      if(eq.hasExplicitDTerm(opId, compId, fieldId))
      {
         computeExplicitTerm<T,SparseMatrix>(eq, opId, compId, rSolverField,  eqStart, fieldId, explicitField, matIdx);
      }
   }

   template <typename T, typename TOperator, typename TData>
      void computeExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx)
   {
      if constexpr((std::is_same<T,MHDFloat>::value || std::is_same<T, MHDComplex>::value ) && (std::is_same<TOperator, SparseMatrixZ>::value && std::is_same<TData, Matrix>::value))
      {
      } else
      {
         // Create pointer to sparse operator
         if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
         {
            ExplicitTermFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(eq, compId, matIdx);
            func.apply(opId, rSolverField, eqStart, fieldId, explicitField);
         }
         else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
         {
            ExplicitTermFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(eq, compId, matIdx);
            func.apply(opId, rSolverField, eqStart, fieldId, explicitField);
         }
         else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
         {
            ExplicitTermFunctor<CouplingIndexType::MODE> func(eq, compId, matIdx);
            func.apply(opId, rSolverField, eqStart, fieldId, explicitField);
         }
         else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
         {
            ExplicitTermFunctor<CouplingIndexType::SINGLE> func(eq, compId, matIdx);
            func.apply(opId, rSolverField, eqStart, fieldId, explicitField);
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EXPLICITTERM_HPP
