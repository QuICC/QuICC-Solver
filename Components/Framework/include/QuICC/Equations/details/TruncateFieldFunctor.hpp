/**
 * @file TruncateFieldFunctor.hpp
 * @brief Implementation of the TruncateField functors
 */

#ifndef QUICC_EQUATIONS_DETAILS_TRUNCATEFIELDFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_TRUNCATEFIELDFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class TruncateFieldFunctor
{
public:
   /**
    * @brief ctor
    */
   TruncateFieldFunctor(const Resolution& res, const Equations::CouplingInformation& cinfo,
      const int matIdx);

   /**
    * @brief deleted default ctor
    */
   TruncateFieldFunctor() = delete;

   /**
    * @brief dtor
    */
   ~TruncateFieldFunctor() = default;

   /**
    * @brief Apply functor
    */
   template <bool IsSet, typename TData, typename TField>
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
    * @brief Matrix index
    */
   const int matIdx;
};

template <CouplingIndexType IndexType>
TruncateFieldFunctor<IndexType>::TruncateFieldFunctor(const Resolution& res, const Equations::CouplingInformation& cinfo,
   const int matIdx) :
    res(res), cinfo(cinfo), matIdx(matIdx)
{
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_TRUNCATEFIELDFUNCTOR_HPP
