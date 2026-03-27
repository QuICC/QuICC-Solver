/**
 * @file StoreSolution.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_STORESOLUTION_HPP
#define QUICC_EQUATIONS_STORESOLUTION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctor.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctorSSR.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctorSMR.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctorM.hpp"
#include "QuICC/Equations/details/StoreSolutionFunctorS.hpp"

namespace QuICC {

namespace Equations {

   template<typename TData, typename TField>
      void storeSolution(TField& field, const Resolution& res, const CouplingInformation& cinfo, const SparseMatrix* pOp, std::shared_ptr<SolutionUpdater> spUp, FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start);

   template<typename TData, typename TField>
      void storeSolution(TField& field, const Resolution& res, const CouplingInformation& cinfo, const SparseMatrix* pOp, std::shared_ptr<SolutionUpdater> spUp, FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      if(cinfo.indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         details::StoreSolutionFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(res, cinfo, pOp, spUp, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(cinfo.indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         details::StoreSolutionFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(res, cinfo, pOp, spUp, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(cinfo.indexType() == CouplingIndexType::MODE)
      {
         details::StoreSolutionFunctor<CouplingIndexType::MODE> func(res, cinfo, pOp, spUp, compId, matIdx);
         func.apply(field, storage, start);
      }
      else if(cinfo.indexType() == CouplingIndexType::SINGLE)
      {
         details::StoreSolutionFunctor<CouplingIndexType::SINGLE> func(res, cinfo, pOp, spUp, compId, matIdx);
         func.apply(field, storage, start);
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_STORESOLUTION_HPP
