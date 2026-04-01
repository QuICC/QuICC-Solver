/**
 * @file AddSource.hpp
 * @brief AddSource implementation
 */

#ifndef QUICC_EQUATIONS_ADDSOURCE_HPP
#define QUICC_EQUATIONS_ADDSOURCE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/details/AddSourceFunctor.hpp"
#include "QuICC/Equations/details/AddSourceFunctorSMR.hpp"
#include "QuICC/Equations/details/AddSourceFunctorSSR.hpp"
#include "QuICC/Equations/details/AddSourceFunctorM.hpp"
#include "QuICC/Equations/details/AddSourceFunctorS.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void addSource(const Resolution& res, const CouplingInformation& cinf0, Spectral::Kernel::SharedISpectralKernel spSrc, const TField& field, TData& storage, const int matIdx, const int start);

   template <typename TData, typename TField> void addSource(const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spSrc, const TField& field, TData& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(cinfo.indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         details::AddSourceFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(res, cinfo, spSrc, matIdx);
         func.apply(field, storage, start);
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(cinfo.indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         details::AddSourceFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(res, cinfo, spSrc, matIdx);
         func.apply(field, storage, start);
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(cinfo.indexType() == CouplingIndexType::MODE)
      {
         details::AddSourceFunctor<CouplingIndexType::MODE> func(res, cinfo, spSrc, matIdx);
         func.apply(field, storage, start);
      }
      // There is a single matrix
      else if(cinfo.indexType() == CouplingIndexType::SINGLE)
      {
         details::AddSourceFunctor<CouplingIndexType::SINGLE> func(res, cinfo, spSrc, matIdx);
         func.apply(field, storage, start);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ADDSOURCE_HPP
