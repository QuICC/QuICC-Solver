/**
 * @file SetBoundaryValueFunctor.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class SetBoundaryValueFunctor
{
public:
   /**
    * @brief ctor
    */
   SetBoundaryValueFunctor(const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spBoundary,
      const int matIdx);

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
   const CouplingInformation& cinfo;

   /**
    * @brief Source kernel
    */
   Spectral::Kernel::SharedISpectralKernel spBoundary;

   /**
    * @brief Matrix index
    */
   const int matIdx;
};

template <CouplingIndexType IndexType>
SetBoundaryValueFunctor<IndexType>::SetBoundaryValueFunctor(
   const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spBoundary,
   const int matIdx) :
    res(res), cinfo(cinfo), spBoundary(spBoundary), matIdx(matIdx)
{}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_SETBOUNDARYVALUEFUNCTOR_HPP
