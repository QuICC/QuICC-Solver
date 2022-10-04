/**
 * @file IEquation.cpp
 * @brief Source of building block for the implementation of a time dependend evolution equation
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/IEquation.hpp"

// Project includes
//
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/PhysicalKernels/DoNothing.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"

namespace QuICC {

namespace Equations {

   IEquation::IEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : EquationData(spEqParams, spScheme, spBackend)
   {
   }

   IEquation::~IEquation()
   {
   }

   void IEquation::init(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Set the coupling
      this->setCoupling();

      // Add the nonlinear integration components
      this->setNLComponents();
   }

   std::vector<Transform::TransformPath> IEquation::forwardPaths()
   {
      std::vector<Transform::TransformPath> paths;

      std::shared_ptr<Transform::ITransformSteps>  spSteps = Transform::createTransformSteps(this->res().sim().spSpatialScheme());

      if(this->requirements(this->name()).isScalar())
      {
         if(this->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear())
         {
            if(this->mForwardPathsType == FWD_IS_FIELD)
            {
               paths = spSteps->forwardScalar(this->nlComponents());

            } else if(this->mForwardPathsType == FWD_IS_NONLINEAR)
            {
               paths = spSteps->forwardNLScalar(this->nlComponents());

            } else
            {
               throw std::logic_error("Custom forward path selected but not defined");
            }
         }
      } else
      {
         if(this->couplingInfo(this->res().sim().ss().spectral().ONE()).hasNonlinear())
         {
            if(this->mForwardPathsType == FWD_IS_FIELD)
            {
               paths = spSteps->forwardVector(this->nlComponents());

            } else if(this->mForwardPathsType == FWD_IS_NONLINEAR)
            {
               paths = spSteps->forwardNLVector(this->nlComponents());

            } else
            {
               throw std::logic_error("Custom forward path selected but not defined");
            }
         }
      }

      return paths;
   }

   void IEquation::initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      //
      // Initialise the galerkin stencils (if activated and required)
      //
      this->initGalerkinStencils(spBcIds, compId);

      //
      // Initialise quasi inverse operator
      //
      this->initQIMatrices(spBcIds, compId);

      //
      // Initialise the explicit linear operators
      //
      this->initExplicitMatrices(spBcIds, compId, ModelOperator::ExplicitLinear::id());

      //
      // Initialise the explicit nonlinear operators
      //
      this->initExplicitMatrices(spBcIds, compId, ModelOperator::ExplicitNonlinear::id());

      //
      // Initialise the explicit nextstep operators
      //
      this->initExplicitMatrices(spBcIds, compId, ModelOperator::ExplicitNextstep::id());
   }

   void IEquation::initGalerkinStencils(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      if(this->couplingInfo(compId).isGalerkin())
      {
         // Get the number of systems
         int nSystems = this->couplingInfo(compId).nSystems();

         this->mGStencils.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator sIt = this->mGStencils.find(compId);
         sIt->second.reserve(nSystems);
         for(int i = 0; i < nSystems; ++i)
         {
            sIt->second.push_back(SparseMatrix());

            this->setGalerkinStencil(compId, sIt->second.back(), i);
         }
      }
   }

   void IEquation::initQIMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      if(this->couplingInfo(compId).hasQuasiInverse())
      {
         // Get the number of systems
         int nSystems = this->couplingInfo(compId).nSystems();

         //
         // Initialise the quasi inverse operators
         //
         SpectralFieldId  fieldId = std::make_pair(this->name(), compId);

         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            // Get block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitBlock(compId, tmpMat.at(i), ModelOperator::ExplicitNonlinear::id(), fieldId, i);

            isComplex = isComplex || (tmpMat.at(i).imag().nonZeros() > 0);
         }

         // Select real or complex operator
         if(isComplex)
         {
            this->mQIZMatrices.insert(std::make_pair(compId, std::vector<SparseMatrixZ>()));
            this->mQIZMatrices.find(compId)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).real().cast<MHDComplex>() + Math::cI*tmpMat.at(i).imag();
               this->mQIZMatrices.find(compId)->second.push_back(tmp);
            }
         } else
         {
            this->mQIDMatrices.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
            this->mQIDMatrices.find(compId)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->mQIDMatrices.find(compId)->second.push_back(SparseMatrix());

               this->mQIDMatrices.find(compId)->second.back().swap(tmpMat.at(i).real());
            }
         }
      }
   }

   void IEquation::initExplicitMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId, const std::size_t opId)
   {
      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

      //
      // Initialise the explicit operators
      //
      CouplingInformation::FieldId_range fRange = this->couplingInfo(compId).explicitRange(opId);
      for(auto fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitBlock(compId, tmpMat.at(i), opId, *fIt, i);

            isComplex = isComplex || (tmpMat.at(i).imag().nonZeros() > 0);
         }

         // Create key
         std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, *fIt);

         // Select real or complex operator
         if(isComplex)
         {
            this->rEZMatrices(opId).insert(std::make_pair(key, std::vector<SparseMatrixZ>()));
            this->rEZMatrices(opId).find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).real().cast<MHDComplex>() + Math::cI*tmpMat.at(i).imag();
               this->rEZMatrices(opId).find(key)->second.push_back(tmp);
            }
         } else
         {
            this->rEDMatrices(opId).insert(std::make_pair(key, std::vector<SparseMatrix>()));
            this->rEDMatrices(opId).find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->rEDMatrices(opId).find(key)->second.push_back(SparseMatrix());

               this->rEDMatrices(opId).find(key)->second.back().swap(tmpMat.at(i).real());
            }
         }
      }
   }

   void IEquation::dispatchCoupling(FieldComponents::Spectral::Id compId, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features, const Resolution& res)
   {
      bool hasNL = features.at(CouplingFeature::Nonlinear);
      bool hasSource = features.at(CouplingFeature::Source);
      bool hasBoundaryValue = features.at(CouplingFeature::BoundaryValue);
      bool allowExplicit = features.at(CouplingFeature::AllowExplicit);

      bool isComplex;
      std::vector<std::pair<std::size_t,FieldComponents::Spectral::Id> >  imFields;
      std::vector<std::pair<std::size_t,FieldComponents::Spectral::Id> >  exLFields;
      std::vector<std::pair<std::size_t,FieldComponents::Spectral::Id> >  exNLFields;
      std::vector<std::pair<std::size_t,FieldComponents::Spectral::Id> >  exNSFields;
      int indexMode;
      auto fId = std::make_pair(this->name(), compId);
      this->backend().equationInfo(isComplex, imFields, exLFields, exNLFields, exNSFields, indexMode, fId, res);

      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(compId,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), compId);

      // Compute effective starting index for local CPU
      int cpuIZero = iZero;
      if(iZero == 1)
      {
         const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
         if(tRes.idx<Dimensions::Data::DAT3D>(0) == 0 && tRes.idx<Dimensions::Data::DAT2D>(0,0) == 0)
         {
            cpuIZero = 1;
         } else
         {
            cpuIZero = 0;
         }
      } else if(iZero > 1)
      {
         throw std::logic_error("Matrix starting index > 1 is not implemented yet!");
      }

      // General setup: equation type? real/complex solver? start from m = ?
      infoIt.first->second.setGeneral(eqType, isComplex, cpuIZero);

      // Set source flag: has source term?
      infoIt.first->second.setSource(hasSource);

      // Set boundary value flag: has boundary value?
      infoIt.first->second.setBoundaryValue(hasBoundaryValue);

      // Set index type: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
      auto idxType = safe_CouplingIndexType_cast(indexMode);
      auto spCoupling = res.sim().ss().createCouplingTools(idxType);
      infoIt.first->second.setIndexType(idxType, spCoupling);

      // Create implicit field coupling
      for(auto fIt = imFields.cbegin(); fIt != imFields.cend(); ++fIt)
      {
         infoIt.first->second.addImplicitField(fIt->first, fIt->second);
      }

      // Create explicit fields
      bool hasQI = false;
      if(allowExplicit)
      {
         // explicit linear
         for(auto fIt = exLFields.cbegin(); fIt != exLFields.cend(); ++fIt)
         {
            infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitLinear::id());
         }

         // explicit nonlinear
         for(auto fIt = exNLFields.cbegin(); fIt != exNLFields.cend(); ++fIt)
         {
            if(!(fIt->first == this->name() && fIt->second == compId))
            {
               infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitNonlinear::id());
            }
         }

         // explicit nextstep
         for(auto fIt = exNSFields.cbegin(); fIt != exNSFields.cend(); ++fIt)
         {
            infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitNextstep::id());
         }

         // Extract quasi inverse
         auto fIt = std::find(exNLFields.begin(), exNLFields.end(), std::make_pair(this->name(), compId));
         if(fIt != exNLFields.end())
         {
            hasQI = true;
         }
      }

      // Set nonlinear flags: has nonlinear term? has quasi-inverse?
      infoIt.first->second.setNonlinear(hasNL, hasNL && hasQI);

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, eqId.second);

      // Get number of matrices
      int nMat = infoIt.first->second.couplingTools().nMat(res);

      // Set field coupling information
      ArrayI tauNs(nMat);
      ArrayI galerkinNs(nMat);
      MatrixI galerkinShifts(nMat, 3);
      ArrayI rhsCols(nMat);
      ArrayI systemNs(nMat);
      this->backend().operatorInfo(tauNs, galerkinNs, galerkinShifts, rhsCols, systemNs, fId, res, infoIt.first->second.couplingTools(), this->bcIds().getTagMap());

      infoIt.first->second.couplingTools().setTauN(tauNs, res);
      infoIt.first->second.couplingTools().setGalerkinN(galerkinNs, res);
      infoIt.first->second.couplingTools().setRhsN(rhsCols, res);
      infoIt.first->second.couplingTools().setSystemN(systemNs, res);
      infoIt.first->second.setSizes(nMat, tauNs, galerkinNs, galerkinShifts, rhsCols, systemNs);
   }

   void  IEquation::dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id compId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs) const
   {
      // Get list of implicit fields
      CouplingInformation::FieldId_range imRange = this->couplingInfo(compId).implicitRange();

      this->backend().modelMatrix(rModelMatrix, opId, imRange, matIdx, bcType, res, eigs, this->bcIds().getTagMap(), this->eqParams().map());
   }

   void IEquation::dispatchGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare) const
   {
      auto fId = std::make_pair(this->name(), compId);
      this->backend().galerkinStencil(mat, fId, matIdx, res, eigs, makeSquare, this->bcIds().getTagMap(), this->eqParams().map());
   }

   void IEquation::dispatchExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs) const
   {
      auto fId = std::make_pair(this->name(), compId);
      this->backend().explicitBlock(mat, fId, opId, fieldId, matIdx, res, eigs, this->bcIds().getTagMap(), this->eqParams().map());
   }

   void IEquation::setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const
   {
      // This implementation should never get called!
      throw std::logic_error("Called dummy implementation of setGalerkinStencil!");
   }

   void IEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const
   {
      // This implementation should never get called!
      throw std::logic_error("Called dummy implementation of setExplicitBlock!");
   }

   Physical::Kernel::SharedIPhysicalKernel IEquation::spNLKernel() const
   {
      assert(this->mspNLKernel);

      return this->mspNLKernel;
   }

   void IEquation::initNLKernel(const bool force)
   {
      if(force || !this->mspNLKernel)
      {
         // Initialize with trivial "do nothing" physical kernel
         this->mspNLKernel = std::make_shared<Physical::Kernel::DoNothing>();
      }
   }

   void  IEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id compId, const int matIdx, const std::size_t bcType) const
   {
      // This implementation should never get called!
      throw std::logic_error("Called dummy implementation of buildModelMatrix!");
   }

   void IEquation::initConstraintKernel()
   {
      for(auto it = this->mConstraintKernel.begin(); it != this->mConstraintKernel.end(); ++it)
      {
         it->second->setResolution(this->spRes());
      }
   }

   void IEquation::setConstraintKernel(FieldComponents::Spectral::Id compId, Spectral::Kernel::SharedISpectralKernel spKernel)
   {
      if(this->mConstraintKernel.count(compId) > 0)
      {
         throw std::logic_error("Source kernel has already been set");
      }

      this->mConstraintKernel.insert(std::make_pair(compId, spKernel));
   }

   Spectral::Kernel::SharedISpectralKernel IEquation::spConstraintKernel(FieldComponents::Spectral::Id compId) const
   {
      assert(this->mConstraintKernel.count(compId) > 0);

      return this->mConstraintKernel.find(compId)->second;
   }

   void IEquation::initSrcKernel()
   {
      for(auto it = this->mSrcKernel.begin(); it != this->mSrcKernel.end(); ++it)
      {
         it->second->setResolution(this->spRes());
      }
   }

   void IEquation::setSrcKernel(FieldComponents::Spectral::Id compId, Spectral::Kernel::SharedISpectralKernel spKernel)
   {
      if(this->mSrcKernel.count(compId) > 0)
      {
         throw std::logic_error("Source kernel has already been set");
      }

      this->mSrcKernel.insert(std::make_pair(compId, spKernel));
   }

   Spectral::Kernel::SharedISpectralKernel IEquation::spSrcKernel(FieldComponents::Spectral::Id compId) const
   {
      assert(this->mSrcKernel.count(compId) > 0);

      return this->mSrcKernel.find(compId)->second;
   }
}
}
