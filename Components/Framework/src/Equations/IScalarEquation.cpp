/**
 * @file IScalarEquation.cpp
 * @brief Source of scalar equation interface
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/Dispatchers.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"
#include "QuICC/Transform/Path/Empty.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"
#include "QuICC/Transform/Path/I2ScalarNl.hpp"

namespace QuICC {

namespace Equations {

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IFieldEquation(spEqParams, spScheme, spBackend)
   {
   }

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions)
      : IFieldEquation(spEqParams, spScheme, spBackend, spOptions)
   {
   }

   void IScalarEquation::setUnknown(Framework::Selector::VariantSharedScalarVariable spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   Framework::Selector::VariantSharedScalarVariable IScalarEquation::spUnknown() const
   {
      // Safety assert
      assert(std::visit([](auto&& p)->bool {return (p != nullptr);}, this->mspUnknown));

      return this->mspUnknown;
   }

   SharedResolution IScalarEquation::spRes() const
   {
      return std::visit([](auto&& p){return p->dom(0).spRes();}, this->spUnknown());
   }

   const Resolution& IScalarEquation::res() const
   {
      return std::visit([](auto&& p)->const Resolution& {return p->dom(0).res();}, this->spUnknown());
   }

   int IScalarEquation::nSpectral() const
   {
      return this->mRequirements.field(this->name()).spectralIds().size();
   }

   typename IScalarEquation::SpectralComponent_range IScalarEquation::spectralRange() const
   {
      return std::make_pair(this->mRequirements.field(this->name()).spectralIds().begin(), this->mRequirements.field(this->name()).spectralIds().end());
   }

   void IScalarEquation::initSpectralMatrices()
   {
      // Make sure it is safe to do nothing
      bool needInit = this->couplingInfo(FieldComponents::Spectral::SCALAR).hasQuasiInverse();

      // Check for Galerkin stencils
      needInit = needInit || this->couplingInfo(FieldComponents::Spectral::SCALAR).isGalerkin();

      // Check explicit linear operators
      CouplingInformation::FieldId_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange(ModelOperator::ExplicitLinear::id());
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);
      // Check explicit nonlinear operators
      fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange(ModelOperator::ExplicitNonlinear::id());
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);
      // Check explicit nextstep operators
      fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange(ModelOperator::ExplicitNextstep::id());
      needInit = needInit || (std::distance(fRange.first,fRange.second) > 0);

      // Initialise spectral matrices
      if(needInit)
      {
         this->initSpectralMatricesComponent(this->mspBcIds, FieldComponents::Spectral::SCALAR);
      }
   }

   void IScalarEquation::defineCoupling(FieldComponents::Spectral::Id compId, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features)
   {
      auto infoIt = this->mCouplingInfos.insert(std::make_pair(compId,CouplingInformation()));
      auto& cinfo = infoIt.first->second;
      auto fieldId = std::make_pair(this->name(), compId);
      dispatchCoupling(cinfo, fieldId, eqType, iZero, features, this->res(), this->backend(), this->bcIds().map());
   }

   void  IScalarEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id compId, const int matIdx, const std::size_t bcType) const
   {
      const auto& cinfo = this->couplingInfo(FieldComponents::Spectral::SCALAR);
      dispatchModelMatrix(rModelMatrix, opId, compId, matIdx, bcType, this->res(), this->backend(), cinfo, this->bcIds().map(), this->eqParams().map());
   }

   void IScalarEquation::setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const
   {
      const auto& cinfo = this->couplingInfo(FieldComponents::Spectral::SCALAR);
      auto fieldId = std::make_pair(this->name(), compId);
      dispatchGalerkinStencil(fieldId, mat, matIdx, this->res(), false, this->backend(), cinfo, this->bcIds().map(), this->eqParams().map());
   }

   void IScalarEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId exId, const int matIdx) const
   {
      const auto& cinfo = this->couplingInfo(FieldComponents::Spectral::SCALAR);
      auto fieldId = std::make_pair(this->name(), compId);
      dispatchExplicitBlock(fieldId, mat, opId, exId, matIdx, this->res(), this->backend(), cinfo, this->bcIds().map(), this->eqParams().map());
   }

   void IScalarEquation::setNLComponents()
   {
      if(this->mForwardPathsType == FWD_IS_NONLINEAR)
      {
         this->addNLComponent(FieldComponents::Spectral::SCALAR, Transform::Path::I2ScalarNl::id());
      }
      else if(this->mForwardPathsType == FWD_IS_FIELD)
      {
         this->addNLComponent(FieldComponents::Spectral::SCALAR, Transform::Path::Scalar::id());
      }
   }

   std::vector<bool> IScalarEquation::disabledBackwardPaths() const
   {
      std::vector<bool> disabled = {false, false, false};

      return disabled;
   }

   std::vector<Transform::TransformPath> IScalarEquation::defaultBackwardPaths(const std::size_t pathId) const
   {
      // Disable some paths
      auto disabled = this->disabledBackwardPaths();
      const bool& disabledPhys = disabled.at(0);
      const bool& disabledGrad = disabled.at(1);
      const bool& disabledGrad2 = disabled.at(2);

      std::vector<Transform::TransformPath> paths;

      auto spSteps = this->transformSteps();

      std::size_t disabledPathId = Transform::Path::Empty::id();

      auto makeMap = [&](auto&& enabled, const bool disabled)
      {
         std::map<typename std::remove_reference<decltype(enabled)>::type::key_type,std::size_t> m;
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
         std::map<FieldComponents::Physical::Id, bool> e = {{FieldComponents::Physical::SCALAR, true}};
         auto compsMap = makeMap(e, disabledPhys);
         auto b = spSteps->backwardScalar(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, this->spUnknown()))
      {
         auto compsMap = std::visit(
               [&](auto&& p)
               {
                  return makeMap(p->dom(0).grad().enabled(), disabledGrad);
               }, this->spUnknown());
         auto b = spSteps->backwardGradient(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, this->spUnknown()))
      {
         auto compsMap = std::visit(
               [&](auto&& p)
               {
                  return makeMap(p->dom(0).grad2().enabled(), disabledGrad2);
               }, this->spUnknown());
         auto b = spSteps->backwardGradient2(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      return paths;
   }

   std::vector<Transform::TransformPath> IScalarEquation::backwardPaths()
   {
      return this->defaultBackwardPaths(Transform::Path::Scalar::id());
   }

   void IScalarEquation::setConstraintKernel(Spectral::Kernel::SharedISpectralKernel spKernel)
   {
      this->setConstraintKernel(FieldComponents::Spectral::SCALAR, spKernel);
   }

   void IScalarEquation::setSrcKernel(Spectral::Kernel::SharedISpectralKernel spKernel)
   {
      this->setSrcKernel(FieldComponents::Spectral::SCALAR, spKernel);
   }

   void IScalarEquation::corruptUnknown(FieldComponents::Spectral::Id compId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      std::visit(
            [&](auto&& p)
            {
               p->rDom(0).rPerturbation().rComp(compId).setConstant(42.42);
            }, this->spUnknown());
   }
} // Equations
} // QuICC
