/**
 * @file IScalarEquation.cpp
 * @brief Source of scalar equation interface
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/IScalarEquation.hpp"

// Project includes
//
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"
#include "QuICC/Transform/Path/I2ScalarNL.hpp"

namespace QuICC {

namespace Equations {

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : IFieldEquation(spEqParams, spScheme, spBackend)
   {
   }

   IScalarEquation::~IScalarEquation()
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

   void IScalarEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features)
   {
      this->dispatchCoupling(comp, eqType, iZero, features, this->res());
   }

   void  IScalarEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType) const
   {
      this->dispatchModelMatrix(rModelMatrix, opId, comp, matIdx, bcType, this->res(), this->couplingInfo(FieldComponents::Spectral::SCALAR).couplingTools().getIndexes(this->res(), matIdx));
   }

   void IScalarEquation::setGalerkinStencil(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      this->dispatchGalerkinStencil(comp, mat, matIdx, this->res(), this->couplingInfo(FieldComponents::Spectral::SCALAR).couplingTools().getIndexes(this->res(), matIdx));
   }

   void IScalarEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const
   {
      this->dispatchExplicitBlock(compId, mat, opId, fieldId, matIdx, this->res(), this->couplingInfo(FieldComponents::Spectral::SCALAR).couplingTools().getIndexes(this->res(), matIdx));
   }

   void IScalarEquation::setNLComponents()
   {
      if(this->mForwardPathsType == FWD_IS_NONLINEAR)
      {
         this->addNLComponent(FieldComponents::Spectral::SCALAR, Transform::Path::I2ScalarNL::id());
      }
      else if(this->mForwardPathsType == FWD_IS_FIELD)
      {
         this->addNLComponent(FieldComponents::Spectral::SCALAR, Transform::Path::Scalar::id());
      }
   }

   std::vector<Transform::TransformPath> IScalarEquation::backwardPaths()
   {
      std::vector<Transform::TransformPath> paths;

      std::shared_ptr<Transform::ITransformSteps>  spSteps = Transform::createTransformSteps(this->res().sim().spSpatialScheme());

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, this->spUnknown()))
      {
         std::map<FieldComponents::Physical::Id,bool> compsMap;
         compsMap.insert(std::make_pair(FieldComponents::Physical::SCALAR, true));
         auto b = spSteps->backwardScalar(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, this->spUnknown()))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<FieldComponents::Physical::Id,bool>{return (p->dom(0).grad().enabled());}, this->spUnknown());
         auto b = spSteps->backwardGradient(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, this->spUnknown()))
      {
         auto compsMap = std::visit([&](auto&& p)->std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>{return (p->dom(0).grad2().enabled());}, this->spUnknown());
         auto b = spSteps->backwardGradient2(compsMap);
         paths.insert(paths.end(), b.begin(), b.end());
      }

      return paths;
   }

   void IScalarEquation::setConstraintKernel(Spectral::Kernel::SharedISpectralKernel spKernel)
   {
      this->setConstraintKernel(FieldComponents::Spectral::SCALAR, spKernel);
   }

   void IScalarEquation::setSrcKernel(Spectral::Kernel::SharedISpectralKernel spKernel)
   {
      this->setSrcKernel(FieldComponents::Spectral::SCALAR, spKernel);
   }
}
}
