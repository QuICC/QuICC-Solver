/**
 * @file IEquation.cpp
 * @brief Source of building block for the implementation of a time dependend evolution equation
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/IEquation.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/PhysicalKernels/DoNothing.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include <unsupported/Eigen/SparseExtra>
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX

namespace QuICC {

namespace Equations {

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
namespace debug {
   /// Create filename to write model operator to MatrixMarket file
   void filenameWriteModelMatrix(const std::string& opName, const std::vector<SpectralFieldId>& tags, const int matIdx);

   /// Write decoupled complex model operator to MatrixMarket file
   void writeModelMatrix(const DecoupledZSparse& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const int matIdx);

   /// Write real model operator to MatrixMarket file
   void writeModelMatrix(const SparseMatrix& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const int matIdx);

   std::string filenameModelMatrix(const std::string& opName, const std::vector<SpectralFieldId>& tags, const int matIdx)
   {
      std::stringstream ss;
      ss << opName;
      for(const auto& f: tags)
      {
         ss << "_" << PhysicalNames::Coordinator::tag(f.first);
         if(f.second == FieldComponents::Spectral::TOR)
         {
            ss <<  "_tor";
         }
         else if(f.second == FieldComponents::Spectral::POL)
         {
            ss <<  "_pol";
         }
      }
      ss << "_" << matIdx;

      return ss.str();
   }

   void writeModelMatrix(const DecoupledZSparse& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const int matIdx)
   {
      auto baseName = filenameModelMatrix(opName, tags, matIdx);

      std::string matName = baseName + "_re.mtx";
      Eigen::saveMarket(mat.real(), matName);
      matName = baseName + "_im.mtx";
      Eigen::saveMarket(mat.imag(), matName);
   }

   void writeModelMatrix(const SparseMatrix& mat, const std::string& opName, const std::vector<SpectralFieldId>& tags, const int matIdx)
   {
      auto baseName = filenameModelMatrix(opName, tags, matIdx);
      std::string matName = baseName + ".mtx";
      Eigen::saveMarket(mat, matName);
   }
}
#endif

   IEquation::IEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : EquationData(spEqParams, spScheme, spBackend)
   {
   }

   IEquation::IEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions)
      : EquationData(spEqParams, spScheme, spBackend, spOptions)
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

      Model::EquationInfo eqInfo;
      auto fId = std::make_pair(this->name(), compId);
      this->backend().equationInfo(eqInfo, fId, res);

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
      infoIt.first->second.setGeneral(eqType, eqInfo.isComplex, cpuIZero, eqInfo.isSplitEquation);

      // Set source flag: has source term?
      infoIt.first->second.setSource(hasSource);

      // Set boundary value flag: has boundary value?
      infoIt.first->second.setBoundaryValue(hasBoundaryValue);

      // Set index type: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
      auto idxType = safe_CouplingIndexType_cast(eqInfo.indexMode);
      auto spCoupling = res.sim().ss().createCouplingTools(idxType);
      infoIt.first->second.setIndexType(idxType, spCoupling);

      // Create implicit field coupling
      for(auto fIt = eqInfo.im.cbegin(); fIt != eqInfo.im.cend(); ++fIt)
      {
         infoIt.first->second.addImplicitField(fIt->first, fIt->second);
      }

      // Create explicit fields
      bool hasQI = false;
      if(allowExplicit)
      {
         // explicit linear
         for(auto fIt = eqInfo.exL.cbegin(); fIt != eqInfo.exL.cend(); ++fIt)
         {
            infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitLinear::id());
         }

         // explicit nonlinear
         for(auto fIt = eqInfo.exNL.cbegin(); fIt != eqInfo.exNL.cend(); ++fIt)
         {
            if(!(fIt->first == this->name() && fIt->second == compId))
            {
               infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitNonlinear::id());
            }
         }

         // explicit nextstep
         for(auto fIt = eqInfo.exNS.cbegin(); fIt != eqInfo.exNS.cend(); ++fIt)
         {
            infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::ExplicitNextstep::id());
         }

         // Extract quasi inverse
         auto fIt = std::find(eqInfo.exNL.begin(), eqInfo.exNL.end(), std::make_pair(this->name(), compId));
         if(fIt != eqInfo.exNL.end())
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
      Model::OperatorInfo opInfo(nMat);
      this->backend().operatorInfo(opInfo, fId, res, infoIt.first->second.couplingTools(), this->bcIds().map());

      infoIt.first->second.couplingTools().setTauN(opInfo.tauN, res);
      infoIt.first->second.couplingTools().setGalerkinN(opInfo.galN, res);
      infoIt.first->second.couplingTools().setRhsN(opInfo.rhsCols, res);
      infoIt.first->second.couplingTools().setSystemN(opInfo.sysN, res);
      infoIt.first->second.setSizes(nMat, opInfo.tauN, opInfo.galN, opInfo.galShift, opInfo.rhsCols, opInfo.sysN);
   }

   void  IEquation::dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id compId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs) const
   {
      // Get list of implicit fields
      CouplingInformation::FieldId_range imRange = this->couplingInfo(compId).implicitRange();

      this->backend().modelMatrix(rModelMatrix, opId, imRange, matIdx, bcType, res, eigs, this->bcIds().map(), this->eqParams().map());

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
      auto opName = ModelOperator::Coordinator::tag(opId);
      auto tags =  std::vector<SpectralFieldId>(imRange.first, imRange.second);
      debug::writeModelMatrix(rModelMatrix, opName, tags, matIdx);
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   }

   void IEquation::dispatchGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare) const
   {
      auto fId = std::make_pair(this->name(), compId);
      this->backend().galerkinStencil(mat, fId, matIdx, res, eigs, makeSquare, this->bcIds().map(), this->eqParams().map());

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
      std::string opName = "galerkin_stencil";
      if(makeSquare)
      {
         opName += "_sq";
      }
      std::vector<SpectralFieldId> tags = {fId};
      debug::writeModelMatrix(mat, opName, tags, matIdx);
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   }

   void IEquation::dispatchExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs) const
   {
      auto fId = std::make_pair(this->name(), compId);
      this->backend().explicitBlock(mat, fId, opId, fieldId, matIdx, res, eigs, this->bcIds().map(), this->eqParams().map());

#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
      auto opName = ModelOperator::Coordinator::tag(opId);
      std::vector<SpectralFieldId> tags = {fId, fieldId};
      debug::writeModelMatrix(mat, opName, tags, matIdx);
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
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

   void IEquation::initConstraintKernel(const std::shared_ptr<std::vector<Array> > spMesh)
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
} // Equations
} // QuICC
