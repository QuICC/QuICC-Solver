/**
 * @file SetBoundaryValueFunctor.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTOR_HPP
#define QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"

namespace QuICC {

namespace Equations {

template <CouplingIndexType IndexType>
class SetBoundaryValueFunctor
{
 public:
   /**
    * @brief ctor
    */
   SetBoundaryValueFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief deleted default ctor
    */
   SetBoundaryValueFunctor() = delete;

   /**
    * @brief dtor
    */
   ~SetBoundaryValueFunctor() = default;

   /**
    * @brief Set boundary value
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void apply(const TField& field, TData& storage, const int start);

 private:
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
SetBoundaryValueFunctor<IndexType>::SetBoundaryValueFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   : eq(&eq), compId(compId), matIdx(matIdx)
{
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETBOUNDARYVALUEFUNCTOR_HPP
