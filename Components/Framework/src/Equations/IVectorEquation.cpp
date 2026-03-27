/**
 * @file IVectorEquation.cpp
 * @brief Source of vector equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Equations/Dispatchers.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"
#include "QuICC/Transform/Path/Empty.hpp"
#include "QuICC/Transform/Path/TorPol.hpp"

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
      auto infoIt = this->mCouplingInfos.insert(std::make_pair(comp,CouplingInformation()));
      auto& cinfo = infoIt.first->second;
      dispatchCoupling(cinfo, this->name(), comp, eqType, iZero, features, this->res(), this->backend(), this->bcIds().map());
   }

   void  IVectorEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType) const
   {
      const auto& cinfo = this->couplingInfo(comp);
      dispatchModelMatrix(rModelMatrix, opId, comp, matIdx, bcType, this->res(), this->backend(), cinfo, this->bcIds().map(), this->eqParams().map());
   }

   void IVectorEquation::setGalerkinStencil(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      const auto& cinfo = this->couplingInfo(comp);
      dispatchGalerkinStencil(this->name(), comp, mat, matIdx, this->res(), false, this->backend(), cinfo, this->bcIds().map(), this->eqParams().map());
   }

   void IVectorEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const
   {
      const auto& cinfo = this->couplingInfo(compId);
      dispatchExplicitBlock(this->name(), compId, mat, opId, fieldId, matIdx, this->res(), this->backend(), cinfo, this->bcIds().map(), this->eqParams().map());
   }

   std::vector<bool> IVectorEquation::disabledBackwardPaths() const
   {
      std::vector<bool> disabled = {false, false, false};

      return disabled;
   }

   std::vector<Transform::TransformPath> IVectorEquation::defaultBackwardPaths(const std::size_t pathId) const
   {
      // Disable some paths
      auto disabled = this->disabledBackwardPaths();
      const bool& disabledPhys = disabled.at(0);
      const bool& disabledGrad = disabled.at(1);
      const bool& disabledCurl = disabled.at(2);

      std::vector<Transform::TransformPath> paths;

      auto spSteps = this->transformSteps();

      const std::size_t disabledPathId = Transform::Path::Empty::id();

      auto makeMap = [&](auto&& enabled, const bool disabled)
      {
         std::map<FieldComponents::Physical::Id,std::size_t> m;
         for(auto&& c: enabled)
         {
            std::size_t id = disabledPathId;
            if(c.second && !disabled)
            {
               id = pathId;
            }
            m.try_emplace(c.first,id);
         }
         return m;
      };

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, this->spUnknown()))
      {
         auto compsMap = std::visit(
               [&](auto&& p)
               {
                  return makeMap(p->dom(0).phys().enabled(), disabledPhys);
               }, this->spUnknown());
         auto branches = spSteps->backwardVector(compsMap);
         paths.insert(paths.end(), branches.begin(), branches.end());
      }

      // vector form
      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad(false));}, this->spUnknown()))
      {
         auto range = this->spectralRange();
         for(auto it = range.first; it != range.second; ++it)
         {
            auto compsMap = std::visit(
               [&](auto&& p)
               {
                  return makeMap(p->dom(0).grad(*it).enabled(), disabledGrad);
               }, this->spUnknown());
            auto b = spSteps->backwardVGradient(*it, compsMap);
            paths.insert(paths.end(), b.begin(), b.end());
         }
      }
      // tensor form
      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad(true));}, this->spUnknown()))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>{return (p->dom(0).grad().enabled());}, this->spUnknown());
         if(disabledGrad)
         {
            for(auto&& c: compsMap)
            {
               c.second = false;
            }
         }
         auto b = spSteps->backwardGradient(compsMap);
         paths.insert(paths.end(),b.begin(), b.end());
      }
      

//      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, this->spUnknown()))
//      {
//          Grad2 is not yet implemented yet implemented
//      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasCurl());}, this->spUnknown()))
      {
         auto compsMap = std::visit(
               [&](auto&& p)
               {
                  return makeMap(p->dom(0).curl().enabled(), disabledCurl);
               }, this->spUnknown());
         auto b = spSteps->backwardCurl(compsMap);
         paths.insert(paths.end(),b.begin(), b.end());
      }

      return paths;
   }

   std::vector<Transform::TransformPath> IVectorEquation::backwardPaths()
   {
      return this->defaultBackwardPaths(Transform::Path::TorPol::id());
   }

   void IVectorEquation::corruptUnknown(FieldComponents::Spectral::Id compId)
   {
      std::visit(
            [&](auto&& p)
            {
               p->rDom(0).rPerturbation().rComp(compId).rData().setConstant(42.42);
            }, this->spUnknown());
   }

} // Equations
} // QuICC
