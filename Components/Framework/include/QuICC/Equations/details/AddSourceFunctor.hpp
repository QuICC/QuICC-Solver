/**
 * @file AddSourceFunctor.hpp
 * @brief AddSourceFunctor implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_ADDSOURCEFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_ADDSOURCEFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class AddSourceFunctor
{
public:
   /**
    * @brief ctor
    */
   AddSourceFunctor(const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spSrc,
      const int matIdx);

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
   const CouplingInformation& cinfo;

   /**
    * @brief Source kernel
    */
   Spectral::Kernel::SharedISpectralKernel spSrc;

   /**
    * @brief Matrix index
    */
   const int matIdx;
};

template <CouplingIndexType IndexType>
AddSourceFunctor<IndexType>::AddSourceFunctor(const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spSrc,
   const int matIdx) :
    res(res), cinfo(cinfo), spSrc(spSrc), matIdx(matIdx)
{}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_ADDSOURCEFUNCTOR_HPP
