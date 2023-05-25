/**
 * @file IVectorEquation.cpp
 * @brief Source of vector equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"

namespace QuICC {

namespace Equations {

   IVectorEquation::IVectorEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IFieldEquation(spEqParams, spScheme, spBackend)
   {
   }

   IVectorEquation::IVectorEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions)
      : IFieldEquation(spEqParams, spScheme, spBackend, spOptions)
   {
   }

   void IVectorEquation::setUnknown(Framework::Selector::VariantSharedVectorVariable spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   Framework::Selector::VariantSharedVectorVariable IVectorEquation::spUnknown() const
   {
      // Safety assert
      assert(std::visit([](auto&& p)->bool {return (p != nullptr);}, this->mspUnknown));

      return this->mspUnknown;
   }

   SharedResolution IVectorEquation::spRes() const
   {
      return std::visit([](auto&& p){return p->dom(0).spRes();}, this->spUnknown());
   }

   const Resolution& IVectorEquation::res() const
   {
      return std::visit([](auto&& p)->const Resolution&{return p->dom(0).res();}, this->spUnknown());
   }

   int IVectorEquation::nSpectral() const
   {
      return this->mRequirements.field(this->name()).spectralIds().size();
   }

   typename IVectorEquation::SpectralComponent_range IVectorEquation::spectralRange() const
   {
      return std::make_pair(this->mRequirements.field(this->name()).spectralIds().begin(), this->mRequirements.field(this->name()).spectralIds().end());
   }

   void IVectorEquation::initSpectralMatrices()
   {
      IVectorEquation::SpectralComponent_range range = this->spectralRange();

      for(auto it = range.first; it != range.second; ++it)
      {
         // Make sure it is safe to do nothing
         bool needInit = this->couplingInfo(*it).hasQuasiInverse();

         // Check for Galerkin stencils
         needInit = needInit || this->couplingInfo(*it).isGalerkin();

         // Check for explicit linear operators
         CouplingInformation::FieldId_range fRange = this->couplingInfo(*it).explicitRange(ModelOperator::ExplicitLinear::id());
         needInit = needInit || (std::distance(fRange.first, fRange.second) > 0);
         // Check for explicit nonlinear operators
         fRange = this->couplingInfo(*it).explicitRange(ModelOperator::ExplicitNonlinear::id());
         needInit = needInit || (std::distance(fRange.first, fRange.second) > 0);
         // Check for explicit nextstep operators
         fRange = this->couplingInfo(*it).explicitRange(ModelOperator::ExplicitNextstep::id());
         needInit = needInit || (std::distance(fRange.first, fRange.second) > 0);

         // Initialise spectral matrices
         if(needInit)
         {
            this->initSpectralMatricesComponent(this->mspBcIds, *it);
         }
      }
   }

   void IVectorEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features)
   {
      this->dispatchCoupling(comp, eqType, iZero, features, this->res());
   }

   void  IVectorEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType) const
   {
      this->dispatchModelMatrix(rModelMatrix, opId, comp, matIdx, bcType, this->res(), this->couplingInfo(comp).couplingTools().getIndexes(this->res(), matIdx));
   }

   void IVectorEquation::setGalerkinStencil(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchGalerkinStencil(comp, mat, matIdx, this->res(), this->couplingInfo(comp).couplingTools().getIndexes(this->res(), matIdx));
   }

   void IVectorEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const
   {
      this->dispatchExplicitBlock(compId, mat, opId, fieldId, matIdx, this->res(), this->couplingInfo(compId).couplingTools().getIndexes(this->res(), matIdx));
   }

   std::vector<Transform::TransformPath> IVectorEquation::backwardPaths()
   {
      std::vector<Transform::TransformPath> paths;

      std::shared_ptr<Transform::ITransformSteps>  spSteps = Transform::createTransformSteps(this->res().sim().spSpatialScheme());

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, this->spUnknown()))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).phys().enabled());}, this->spUnknown());
         auto branches = spSteps->backwardVector(compsMap);
         paths.insert(paths.end(), branches.begin(), branches.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, this->spUnknown()))
      {
         auto range = this->spectralRange();
         for(auto it = range.first; it != range.second; ++it)
         {
            auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).grad(*it).enabled());}, this->spUnknown());
            auto b = spSteps->backwardVGradient(*it, compsMap);
            paths.insert(paths.end(), b.begin(), b.end());
         }
      }

// Not yet implemented
//      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, this->spUnknown()))
//      {
//         auto range = this->spectralRange();
//         for(auto it = range.first; it != range.second; ++it)
//         {
//            auto compsMap = std::visit([&](auto&& p)->std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>{return (p->dom(0).grad2(*it).enabled());}, this->spUnknown());
//            auto b = spSteps->backwardVGradient2(*it, compsMap);
//            paths.insert(paths.end(),b.begin(), b.end());
//         }
//      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasCurl());}, this->spUnknown()))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).curl().enabled());}, this->spUnknown());
         auto b = spSteps->backwardCurl(compsMap);
         paths.insert(paths.end(),b.begin(), b.end());
      }

      return paths;
   }

} // Equations
} // QuICC
