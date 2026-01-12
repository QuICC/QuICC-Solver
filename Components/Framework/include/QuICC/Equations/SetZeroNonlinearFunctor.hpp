/**
 * @file SetZeroNonlinearFunctor.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_SETZERONONLINEARFUNCTOR_HPP
#define QUICC_EQUATIONS_SETZERONONLINEARFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"

namespace QuICC {

namespace Equations {

template <CouplingIndexType IndexType>
class SetZeroNonlinearFunctor
{
 public:
   /**
    * @brief ctor
    */
   SetZeroNonlinearFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief deleted default ctor
    */
   SetZeroNonlinearFunctor() = delete;

   /**
    * @brief dtor
    */
   ~SetZeroNonlinearFunctor() = default;

   /**
    * @brief Set nonlinear spectral values to zero
    *
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
SetZeroNonlinearFunctor<IndexType>::SetZeroNonlinearFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   : eq(&eq), compId(compId), matIdx(matIdx)
{
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETZERONONLINEARFUNCTOR_HPP
