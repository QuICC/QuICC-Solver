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
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#endif
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "Types/Math.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/PhysicalKernels/DoNothing.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"

namespace QuICC {

namespace Equations {

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

      // Set solution updater
      this->initSolutionUpdater();
   }

   std::shared_ptr<Transform::ITransformSteps> IEquation::transformSteps() const
   {
      auto  spSteps = Transform::createTransformSteps(this->res().sim().spSpatialScheme());

      return spSteps;
   }

   std::vector<Transform::TransformPath> IEquation::forwardPaths()
   {
      std::vector<Transform::TransformPath> paths;

      auto spSteps = this->transformSteps();

      if(this->requirements(this->name()).isScalar())
      {
         if(this->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear())
         {
            if(this->mForwardPathsType == FWD_IS_FIELD)
            {
               paths = spSteps->forwardScalar(this->nlComponents());

            }
            else if(this->mForwardPathsType == FWD_IS_NONLINEAR)
            {
               paths = spSteps->forwardNLScalar(this->nlComponents());

            }
            else
            {
               throw std::logic_error("Custom forward path selected but not defined");
            }
         }
      }
      else
      {
         if(this->couplingInfo(this->res().sim().ss().spectral().ONE()).hasNonlinear())
         {
            if(this->mForwardPathsType == FWD_IS_FIELD)
            {
               paths = spSteps->forwardVector(this->nlComponents());

            }
            else if(this->mForwardPathsType == FWD_IS_NONLINEAR)
            {
               paths = spSteps->forwardNLVector(this->nlComponents());

            }
            else
            {
               throw std::logic_error("Custom forward path selected but not defined");
            }
         }
      }

      return paths;
   }

   void IEquation::initSolutionUpdater()
   {
      auto range = this->spectralRange();

      for(auto it = range.first; it != range.second; ++it)
      {
         auto spUp = std::make_shared<SolutionUpdater>();
         this->mSolUps.emplace(*it, spUp);
      }
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
      if(this->mConstraintKernel.count(compId) > 0)
      {
         return this->mConstraintKernel.at(compId);
      }
      else
      {
         return nullptr;
      }
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

   void IEquation::writeDiagnostics(const bool, const bool) const
   {
   }

   void IEquation::linkEquation(std::shared_ptr<IEquation> spEq)
   {
      // Default does nothing
   }

   void IEquation::updateConstraintKernel(const MHDFloat time, const MHDFloat timestep, const bool isFinished)
   {
      DebuggerMacro_msg("(Nothing to update for " + PhysicalNames::Coordinator::tag(this->name()) + "(it = " + std::to_string(this->options().it()) + "))", 5);
      // Default does nothing
   }
} // Equations
} // QuICC
