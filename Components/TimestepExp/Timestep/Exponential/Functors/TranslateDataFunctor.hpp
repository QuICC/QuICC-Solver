/**
 * @file TranslateDataFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "QuICC/IteratorRange.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Translate equation data 
 */
class TranslateDataFunctor
{
   public:
      /**
       * @brief ctor
       */
      TranslateDataFunctor(std::shared_ptr<FunctorData> spData) : mspData(spData) {};

      /**
       * @brief dtor
       */
      ~TranslateDataFunctor() = default;

      template <typename TEqIt>
      void operator()(const SpectralFieldId& id, TEqIt& eqIt);

   private:
      /**
       * @brief Data storage
       */
      std::shared_ptr<FunctorData> mspData;
};

template <typename TEqIt>
void TranslateDataFunctor::operator()(const SpectralFieldId& myId, TEqIt& eqIt)
{
   assert(this->mspData);
   auto& data = *this->mspData;
   auto& ciMap = data.cInfos;

   // Set Resolution
   if(data.spRes== nullptr)
   {
      data.spRes = eqIt->spRes();
   }

   // Set BcIdMap
   if(data.bcIdMap.size() == 0)
   {
      data.bcIdMap = eqIt->bcIds().map();
   }

   // Set EqParamsMap
   if(data.eqParamsMap.size() == 0)
   {
      data.eqParamsMap = eqIt->eqParams().map();
   }

   // Set backend
   if(data.spBackend == nullptr)
   {
      data.spBackend = eqIt->spBackend();
   }

   // Equation information
   data.eqInfos.emplace_back(eqIt->options().it(), eqIt->solveTiming(), myId);

   // Get coupling information
   ciMap.emplace(myId, eqIt->couplingInfo(myId.second));
   // Get field component pointer
   std::visit(
         [&](auto&& p)
         {
            auto&& f = p->rDom(0).rPerturbation().rComp(myId.second);
            using T = std::decay_t<decltype(f)>;
            if constexpr(std::is_same_v<T, Framework::Selector::ComplexScalarField>)
            {
               data.fields.emplace(myId, &f);
            }
            else
            {
               throw std::logic_error("Real scalar fields are not implemented");
            }
   }, eqIt->spUnknown());

   const auto& cinfo = ciMap.at(myId);

   // Set stencils
   if(cinfo.isGalerkin())
   {
      data.stencils.emplace(myId, eqIt->galerkinStencils(myId.second));
   }

   // Set Constraint
   auto spConstraint = eqIt->spConstraintKernel(myId.second);
   if(spConstraint)
   {
      data.constraints.emplace(myId, eqIt->spConstraintKernel(myId.second));
   }

   // Set Boundary value
   if(cinfo.hasBoundaryValue())
   {
      data.bcvalues.emplace(myId, eqIt->spBoundaryKernel(myId.second));
   }

   // Set Source
   if(cinfo.hasSource())
   {
      data.sources.emplace(myId, eqIt->spSourceKernel(myId.second));
   }

   // Set Solution updater
   data.solups.emplace(myId, eqIt->spSolutionUpdater(myId.second));
   
   // Get explicit operators
   std::size_t opId = ModelOperator::ExplicitLinear::id();
   auto r = make_range(cinfo.explicitRange(opId));
   std::map<SpectralFieldId, std::vector<SparseMatrix>> tmpExD;
   std::map<SpectralFieldId, std::vector<SparseMatrixZ>> tmpExZ;
   for(auto& exId: r)
   {
      // Complex linear operator
      if (eqIt->hasExplicitZTerm(opId, myId.second, exId))
      {
         tmpExZ.emplace(exId, eqIt->template explicitOperators<SparseMatrixZ>(opId,
               myId.second, exId));
      }

      // Real linear operator
      if (eqIt->hasExplicitDTerm(opId, myId.second, exId))
      {
         tmpExD.emplace(exId, eqIt->template explicitOperators<SparseMatrix>(opId,
               myId.second, exId));
      }
   }
   if(tmpExZ.size() > 0)
   {
      auto&& exTerm = data.exZTerm;
      if(exTerm.count(opId) == 0)
      {
         exTerm.emplace(opId, std::map<SpectralFieldId, std::map<SpectralFieldId, std::vector<SparseMatrixZ>>>());
      }
      //exTerm.at(opId).emplace(myId, tmpExZ);
   }
   if(tmpExD.size() > 0)
   {
      auto&& exTerm = data.exDTerm;
      if(exTerm.count(opId) == 0)
      {
         exTerm.emplace(opId, std::map<SpectralFieldId, std::map<SpectralFieldId, std::vector<SparseMatrix>>>());
      }
      exTerm.at(opId).emplace(myId, tmpExD);
   }
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
