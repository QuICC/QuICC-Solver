/**
 * @file SetZeroNonlinearFunctor.hpp
 * @brief SetZeroNonlinear implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class SetZeroNonlinearFunctor
{
public:
   /**
    * @brief ctor
    */
   SetZeroNonlinearFunctor(const Resolution& res, const CouplingInformation& cinfo,
      FieldComponents::Spectral::Id compId, const int matIdx);

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
   template <typename TData, typename TField>
   void apply(const TField& field, TData& storage, const int start);

private:
   /**
    * @brief Resolution
    */
   const Resolution& res;

   /**
    * @brief Coupling information
    */
   const Equations::CouplingInformation& cinfo;

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
SetZeroNonlinearFunctor<IndexType>::SetZeroNonlinearFunctor(
   const Resolution& res, const CouplingInformation& cinfo, FieldComponents::Spectral::Id compId,
   const int matIdx) :
    res(res), cinfo(cinfo), compId(compId), matIdx(matIdx)
{}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETZERONONLINEARFUNCTOR_HPP
