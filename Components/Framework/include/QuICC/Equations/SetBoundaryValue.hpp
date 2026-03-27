/**
 * @file SetBoundaryValue.hpp
 * @brief SetBoundaryValue implementation
 */

#ifndef QUICC_EQUATIONS_SETBOUNDARYVALUE_HPP
#define QUICC_EQUATIONS_SETBOUNDARYVALUE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctor.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctorSSR.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctorSMR.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctorM.hpp"
#include "QuICC/Equations/details/SetBoundaryValueFunctorS.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Set boundary value
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void setBoundaryValue(const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spBoundary, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TData, typename TField> void setBoundaryValue(const Resolution& res, const CouplingInformation& cinfo, Spectral::Kernel::SharedISpectralKernel spBoundary, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Set boundary value if required
      if(cinfo.hasBoundaryValue())
      {
         if(cinfo.isGalerkin())
         {
            throw std::logic_error("Galerkin expansion cannot have a nonzero boundary value!");
         }

         // matIdx is the index of the slowest varying direction with a single RHS
         if(cinfo.indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
         {
            details::SetBoundaryValueFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(res, cinfo, spBoundary, compId, matIdx);
            func.apply(field, storage, start);
         }
         // matIdx is the index of the slowest varying direction with multiple RHS
         else if(cinfo.indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
         {
            details::SetBoundaryValueFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(res, cinfo, spBoundary, compId, matIdx);
            func.apply(field, storage, start);
         }
         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         else if(cinfo.indexType() == CouplingIndexType::MODE)
         {
            details::SetBoundaryValueFunctor<CouplingIndexType::MODE> func(res, cinfo, spBoundary, compId, matIdx);
            func.apply(field, storage, start);
         }
         // There is a single matrix
         else if(cinfo.indexType() == CouplingIndexType::SINGLE)
         {
            details::SetBoundaryValueFunctor<CouplingIndexType::SINGLE> func(res, cinfo, spBoundary, compId, matIdx);
            func.apply(field, storage, start);
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_SETBOUNDARYVALUE_HPP
