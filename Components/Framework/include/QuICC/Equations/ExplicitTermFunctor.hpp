/**
 * @file ExplicitTermFunctor.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_EXPLICITTERMFUNCTOR_HPP
#define QUICC_EQUATIONS_EXPLICITTERMFUNCTOR_HPP

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

namespace QuICC {

namespace Equations {

template <CouplingIndexType IndexType>
class ExplicitTermFunctor
{
 public:
   /**
    * @brief ctor
    */
   ExplicitTermFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief deleted default ctor
    */
   ExplicitTermFunctor() = delete;

   /**
    * @brief dtor
    */
   ~ExplicitTermFunctor() = default;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param fieldId       Physical field ID
    * @param explicitField Explicit linear field values
    */
   template <typename T, typename TData>
      void apply(const std::size_t opId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField);

 private:
   template <typename T, typename TOperator,typename TData> void compute(const std::size_t opId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Framework::Selector::ScalarField<T>& explicitField);

   /**
    * @brief Reference to equation
    */
   const IFieldEquation* eq;

   /**
    * @brief Field component ID
    */
   FieldComponents::Spectral::Id compId;

   /**
    * @brief Matrix index
    */
   const int matIdx;
};

   template <CouplingIndexType IndexType>
ExplicitTermFunctor<IndexType>::ExplicitTermFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   : eq(&eq), compId(compId), matIdx(matIdx)
{
}

template <CouplingIndexType IndexType>
   template <typename T, typename TData>
void ExplicitTermFunctor<IndexType>::apply(const std::size_t opId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField)
{
   // Compute with complex linear operator
   if(eq->hasExplicitZTerm(opId, compId, fieldId))
   {
      compute<T,SparseMatrixZ>(opId, rSolverField, eqStart, fieldId, explicitField);
   }

   // Compute with real linear operator
   if(eq->hasExplicitDTerm(opId, compId, fieldId))
   {
      compute<T,SparseMatrix>(opId, rSolverField, eqStart, fieldId, explicitField);
   }
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EXPLICITTERMFUNCTOR_HPP
