/**
 * @file AddSourceFunctor.hpp
 * @brief AddSourceFunctor implementation
 */

#ifndef QUICC_EQUATIONS_ADDSOURCEFUNCTOR_HPP
#define QUICC_EQUATIONS_ADDSOURCEFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"

namespace QuICC {

namespace Equations {

template <CouplingIndexType IndexType>
class AddSourceFunctor
{
 public:
   /**
    * @brief ctor
    */
   AddSourceFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief deleted default ctor
    */
   AddSourceFunctor() = delete;

   /**
    * @brief dtor
    */
   ~AddSourceFunctor() = default;

   /**
    * @brief Add source term
    *
    * @param field   Field
    * @param storage Storage for the equation values
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
AddSourceFunctor<IndexType>::AddSourceFunctor(const IFieldEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   : eq(&eq), compId(compId), matIdx(matIdx)
{
}

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCEFUNCTOR_HPP
