/**
 * @file Coordinator.cpp
 * @brief Source of the high level pseudospectral coordinator
 */

// System includes
//
#include <algorithm>
#include <stdexcept>

// Project includes
//
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/SolveTiming/Before.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Variables/RequirementTools.hpp"
#include "QuICC/TransformCoordinators/TransformCoordinatorTools.hpp"
#include "QuICC/Equations/Tools/EquationTools.hpp"
#include "QuICC/Simulation/SimulationIoTools.hpp"
#include "QuICC/PseudospectralTag/Trivial.hpp"
#include "QuICC/PseudospectralTag/Diagnostic.hpp"
#include "QuICC/PseudospectralTag/Prognostic.hpp"
#include "QuICC/PseudospectralTag/Uninitialized.hpp"
#include "QuICC/PseudospectralTag/Wrapper.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Pseudospectral {

   const std::set<int>& Coordinator::it() const
   {
      return this->mIt;
   }

   MHDFloat Coordinator::time() const
   {
      return this->mTimestepCoordinator.time();
   }

   MHDFloat Coordinator::timestep() const
   {
      return this->mTimestepCoordinator.timestep();
   }

   MHDFloat Coordinator::startTime() const
   {
      return this->mDiagnostics.startTime();
   }

   MHDFloat Coordinator::startTimestep() const
   {
      return this->mDiagnostics.startTimestep();
   }

   void Coordinator::printInfo(std::ostream& stream)
   {
      this->mTimestepCoordinator.printInfo(stream);
   }

   void Coordinator::addEquation(Equations::SharedIScalarEquation spEq, const int it)
   {
      this->addEquation(spEq, PseudospectralTag::Uninitialized::id(), it);
   }

   void Coordinator::addEquation(Equations::SharedIScalarEquation spEq, const std::size_t eqId, const int it)
   {
      DebuggerMacro_showValue("Adding scalar equation for " + PhysicalNames::Coordinator::tag(spEq->name()) +" to pseudospectral coordinator iteration ", 1, it);

      // Set equation parameters if not yet set
      if(!this->mspEqParams)
      {
         this->mspEqParams = spEq->spEqParams();
      }

      if(this->mScalarEquations.count(it) == 0)
      {
         this->mScalarEquations.insert(std::make_pair(it, std::vector<SharedIScalarEquation>()));
      }
      this->mScalarEquations.at(it).push_back(spEq);

      this->mIt.insert(it);

      std::pair<std::size_t,int> key = std::make_pair(eqId, it);
      if(this->mScalarEqMap.count(key) == 0)
      {
         this->mScalarEqMap.insert(std::pair(key, std::vector<SharedIScalarEquation>()));
      }

      this->mScalarEqMap.at(key).push_back(spEq);
      DebuggerMacro_msg("... done", 1);
   }

   void Coordinator::addEquation(Equations::SharedIVectorEquation spEq, const int it)
   {
      this->addEquation(spEq, PseudospectralTag::Uninitialized::id(), it);
   }

   void Coordinator::addEquation(Equations::SharedIVectorEquation spEq, const std::size_t eqId, const int it)
   {
      DebuggerMacro_showValue("Adding vector equation for " + PhysicalNames::Coordinator::tag(spEq->name()) +" to pseudospectral coordinator iteration ", 1, it);

      // Set equation parameters if not yet set
      if(!this->mspEqParams)
      {
         this->mspEqParams = spEq->spEqParams();
      }

      if(this->mVectorEquations.count(it) == 0)
      {
         this->mVectorEquations.insert(std::make_pair(it, std::vector<SharedIVectorEquation>()));
      }
      this->mVectorEquations.at(it).push_back(spEq);

      this->mIt.insert(it);

      std::pair<std::size_t,int> key = std::make_pair(eqId, it);
      if(this->mVectorEqMap.count(key) == 0)
      {
         this->mVectorEqMap.insert(std::pair(key, std::vector<SharedIVectorEquation>()));
      }

      this->mVectorEqMap.at(key).push_back(spEq);
      DebuggerMacro_msg("... done", 1);
   }

   Coordinator::ScalarEquation_range Coordinator::scalarRange(const std::size_t eqId, const int it)
   {
      auto key = std::make_pair(eqId, it);
      assert(this->mScalarEqMap.count(key) > 0);
      auto eqIt = this->mScalarEqMap.find(key);
      auto r = std::make_pair(eqIt->second.begin(), eqIt->second.end());

      return r;
   }

   Coordinator::VectorEquation_range Coordinator::vectorRange(const std::size_t eqId, const int it)
   {
      auto key = std::make_pair(eqId, it);
      assert(this->mVectorEqMap.count(key) > 0);
      auto eqIt = this->mVectorEqMap.find(key);
      auto r = std::make_pair(eqIt->second.begin(), eqIt->second.end());

      return r;
   }

   Coordinator::SharedIScalarEquation Coordinator::scalarEq(const std::size_t eqId, const int it, const int i)
   {
      auto key = std::make_pair(eqId, it);
      assert(this->mScalarEqMap.count(key) > 0);
      auto eqIt = this->mScalarEqMap.find(key);

      assert(eqIt->second.size() > 0);
      return eqIt->second.at(i);
   }

   Coordinator::SharedIVectorEquation Coordinator::vectorEq(const std::size_t eqId, const int it, const int i)
   {
      auto key = std::make_pair(eqId, it);
      assert(this->mVectorEqMap.count(key) > 0);
      auto eqIt = this->mVectorEqMap.find(key);

      assert(eqIt->second.size() > 0);
      return eqIt->second.at(i);
   }

   void Coordinator::evolve()
   {
      DebuggerMacro_msg("Evolving equations", 1);
      Profiler::RegionFixture fix("evolve");

      bool isIntegrating = true;
      while(isIntegrating)
      {
         DebuggerMacro_msg("Time integration sub-step", 2);
         // Update equation time
         this->updateEquationTime(this->time(), this->mTimestepCoordinator.finishedStep());

         // Loop over sub-steps
         for(auto j: this->it())
         {
            DebuggerMacro_showValue("Equation sub-iteration ", 3, j);

            // Compute explicit linear terms
            DebuggerMacro_msg("Explicit equations", 4);
            this->explicitEquations(j);

            QuICCTimer().stop();
            QuICCTimer().update(ExecutionTimer::RUN);
            QuICCTimer().start();

            // Compute the nonlinear terms
            DebuggerMacro_msg("Transform loop", 4);
            this->computeNonlinear(j);

            QuICCTimer().stop();
            QuICCTimer().update(ExecutionTimer::NONLINEAR);
            QuICCTimer().update(ExecutionTimer::RUN);
            QuICCTimer().start();

            // Timestep the equations
            DebuggerMacro_msg("Solve equations", 4);
            this->solveEquations(j);

            QuICCTimer().stop();
            QuICCTimer().update(ExecutionTimer::TIMESTEP);
            QuICCTimer().update(ExecutionTimer::RUN);
            QuICCTimer().start();
         }

         // Update the equations
         for(auto j: this->it())
         {
            DebuggerMacro_msg("Update equations", 4);
            this->updateEquations(j, this->mTimestepCoordinator.finishedStep());
         }

         // Finish timestep iteration
         this->finalizeTimestep();

         isIntegrating = !this->mTimestepCoordinator.finishedStep();
      }
   }

   void Coordinator::initTransformCoordinator()
   {
      // Extract the run options for the equation parameters
      std::map<std::size_t,NonDimensional::SharedINumber> runOptions(this->mspEqParams->map());

      // Compute pack sizes
      std::vector<ArrayI> packs;
      Transform::TransformCoordinatorTools::computePacks(packs, this->mspFwdGrouper, this->mspBwdGrouper, this->mFwdTree, mBwdTree, this->it(), this->mspRes);

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(this->mTransformCoordinator, this->mspFwdGrouper, this->mspBwdGrouper, packs, this->mspRes, runOptions);
   }

   void Coordinator::initParallel(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      StageTimer stage;
      stage.start("intializing transform grouper");

      // Store the shared resolution object
      this->mspRes = spRes;

      // Initialise the transform grouper
      Parallel::setGrouper(descr, this->mspFwdGrouper, this->mspBwdGrouper);

      // Initialise the transform grouper for imposed fields
      Parallel::setGrouper(descr, this->mspImposedFwdGrouper, this->mspImposedBwdGrouper);

      stage.done();
   }

   void Coordinator::updatePhysical(const int it)
   {
      Profiler::RegionFixture fix("updatePhysical");

      // Define transform trees
      this->mTransformCoordinator.defineBwdTransforms(this->mBwdTree.at(it));

      // Compute physical values
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
   }

   void Coordinator::updateSpectral(const int it)
   {
      Profiler::RegionFixture fix("updateSpectral");

      // Define transform trees
      this->mTransformCoordinator.defineFwdTransforms(this->mFwdTree.at(it));

      this->mspFwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mPhysicalKernels.at(it), this->mTransformCoordinator);
   }

   void Coordinator::updateSpectral(const bool isTrivial, const bool isDiagnostic, const bool isPrognostic, const bool isWrapper, const int it)
   {
      std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel> kernels;

      // Get kernels from trivial equations
      if(isTrivial)
      {
         auto sR = this->scalarRange(PseudospectralTag::Trivial::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Trivial::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Get kernels from diagnostic equations
      if(isDiagnostic)
      {
         auto sR = this->scalarRange(PseudospectralTag::Diagnostic::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Diagnostic::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Get kernels from prognostic equations
      if(isPrognostic)
      {
         auto sR = this->scalarRange(PseudospectralTag::Prognostic::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Prognostic::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Get kernels from wrapper equations
      if(isWrapper)
      {
         auto sR = this->scalarRange(PseudospectralTag::Wrapper::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Wrapper::id(), it);
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Set mesh for kernels
      Equations::Tools::setupPhysicalKernels(kernels, this->mTransformCoordinator.mesh());

      //this->mspFwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, kernels, this->mTransformCoordinator);
   }

   void Coordinator::init(const Array& tstep, const SharedSimulationBoundary spBcs)
   {
      StageTimer stage;

      // Make sure types are initialized for both equation types
      for(auto j: this->it())
      {
         // Create empty scalar vector for iteration
         if(this->mScalarEquations.count(j) == 0)
         {
            this->mScalarEquations.insert(std::make_pair(j, std::vector<SharedIScalarEquation>()));
         }

         // Create empty vector vector for iteration
         if(this->mVectorEquations.count(j) == 0)
         {
            this->mVectorEquations.insert(std::make_pair(j, std::vector<SharedIVectorEquation>()));
         }
      }

      stage.start("initializing variables");

      // Initialise the variables and set general variable requirements
      VariableRequirement varInfo;
      for(auto j: this->it())
      {
         assert(this->mScalarEquations.count(j) == 1);
         assert(this->mVectorEquations.count(j) == 1);

         RequirementTools::mergeRequirements(varInfo, this->mScalarEquations.at(j), this->mVectorEquations.at(j));
      }
      RequirementTools::initVariables(this->mScalarVariables, this->mVectorVariables, varInfo, this->mspRes);

      // Map variables to the equations and set nonlinear requirements
      std::vector<std::size_t> unmapped;
      for(auto j: this->it())
      {
         RequirementTools::mapEquationVariables(this->mScalarEquations.at(j), this->mVectorEquations.at(j), this->mScalarVariables, this->mVectorVariables, spBcs, unmapped);
      }

      if(unmapped.size() > 0)
      {
         throw std::logic_error("Variables not mapped to equations are currently not implemented");
      }

      // Initialize Imposed fields
      this->initImposed();

      // Transform trees
      for(auto j: this->it())
      {
         // Build backward tree
         this->mBwdTree.insert(std::make_pair(j, std::vector<Transform::TransformTree>()));
         RequirementTools::buildBackwardTree(this->mBwdTree.at(j), this->mScalarEquations.at(j), this->mVectorEquations.at(j));

         // Build forward tree
         this->mFwdTree.insert(std::make_pair(j, std::vector<Transform::TransformTree>()));
         RequirementTools::buildForwardTree(this->mFwdTree.at(j), this->mScalarEquations.at(j), this->mVectorEquations.at(j));

         // Get nonlinear kernels
         this->mPhysicalKernels.insert(std::make_pair(j, std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel>()));
         // ... for scalar equations
         Equations::Tools::getNonlinearKernels(this->mPhysicalKernels.at(j), this->mScalarEquations.at(j).begin(), this->mScalarEquations.at(j).end());
         // ... for vector equations
         Equations::Tools::getNonlinearKernels(this->mPhysicalKernels.at(j), this->mVectorEquations.at(j).begin(), this->mVectorEquations.at(j).end());
      }

      stage.done();

      // Initialise the transform coordinator
      this->initTransformCoordinator();

      // Setup nonlinear kernels
      for(auto j: this->it())
      {
         Equations::Tools::setupPhysicalKernels(this->mPhysicalKernels.at(j), this->mTransformCoordinator.mesh());
      }

      stage.start("setup equations");

      // Initialise the equations (generate operators, etc)
      this->setupEquations();

      // Sort the equations by type: time/solver/trivial
      this->sortEquations();

      stage.done();
      stage.start("initializing diagnostics");

      // Initialise the diagnostics
      this->mDiagnostics.init(this->transformCoordinator().mesh(), this->scalarVariables(), this->vectorVariables(), tstep, this->mspEqParams->map());
      stage.done();
   }

   void Coordinator::useStateTime(const MHDFloat time, const MHDFloat timestep)
   {
      this->mDiagnostics.useStateTime(time, timestep);
   }

   void Coordinator::prepareEvolution(const std::size_t schemeId)
   {
      // Update equation time
      this->updateEquationTime(this->mDiagnostics.startTime(), false);

      // Initialise all values (solve and nonlinear computations except timestep)
      this->preSolveEquations();

      // Update CFL condition
      this->mDiagnostics.initialCfl();

      // Synchronise diagnostics
      this->mDiagnostics.synchronize();

      // Init timestepper using clf/100 as starting timestep
      std::vector<SharedIScalarEquation> sEqs;
      std::vector<SharedIVectorEquation> vEqs;
      for(auto j: this->it())
      {
         auto sP = this->scalarRange(PseudospectralTag::Prognostic::id(), j);
         sEqs.insert(sEqs.end(), sP.first, sP.second);
         auto vP = this->vectorRange(PseudospectralTag::Prognostic::id(), j);
         vEqs.insert(vEqs.end(), vP.first, vP.second);
      }
      auto sP = std::make_pair(sEqs.begin(), sEqs.end());
      auto vP = std::make_pair(vEqs.begin(), vEqs.end());
      this->mTimestepCoordinator.init(schemeId, this->mDiagnostics.startTime(), this->mDiagnostics.cfl(), this->mDiagnostics.maxError(), sP, vP);

      // Compute physical space values if required
      this->mspImposedBwdGrouper->transform(this->mImposedScalarVariables, this->mImposedVectorVariables, *this->mspImposedTransformCoordinator);
   }

   void Coordinator::cleanupForRun()
   {
      StageTimer stage;
      stage.start("Cleanup equation backends");

      // Loop over iterations
      for(auto j: this->it())
      {
         // Loop over all scalar equations
         auto sit = this->mScalarEquations.find(j);
         for(auto scalEqIt = sit->second.begin(); scalEqIt < sit->second.end(); ++scalEqIt)
         {
            (*scalEqIt)->cleanupBackend();
         }

         // Loop over all vector equations
         auto vit = this->mVectorEquations.find(j);
         for(auto vectEqIt = vit->second.begin(); vectEqIt < vit->second.end(); ++vectEqIt)
         {
            (*vectEqIt)->cleanupBackend();
         }
      }

      // Cleanup imposed transform grouper
      this->mspImposedFwdGrouper.reset();
      this->mspImposedBwdGrouper.reset();
      this->mspImposedTransformCoordinator.reset();

      stage.done();
   }

   void Coordinator::initImposed()
   {
      StageTimer stage;

      stage.start("initializing imposed variables");

      // Initialise the variables and set general variable requirements
      VariableRequirement varInfo;
      for(auto j: this->it())
      {
         assert(this->mScalarEquations.count(j) == 1);
         assert(this->mVectorEquations.count(j) == 1);

         RequirementTools::mergeImposedRequirements(varInfo, this->mScalarEquations.at(j), this->mVectorEquations.at(j));
      }
      RequirementTools::initVariables(this->mImposedScalarVariables, this->mImposedVectorVariables, varInfo, this->mspRes);

      // Map variables to the equations and set nonlinear requirements
      std::vector<std::size_t> unmapped;
      for(auto j: this->it())
      {
         RequirementTools::mapImposedVariables(this->mScalarEquations.at(j), this->mVectorEquations.at(j), this->mImposedScalarVariables, this->mImposedVectorVariables, unmapped);
      }

      if(unmapped.size() > 0)
      {
         throw std::logic_error("Variables not mapped to equations are currently not implemented");
      }

      // Transform trees
      std::vector<Transform::TransformTree> forwardTree;
      std::vector<Transform::TransformTree> backwardTree;
      RequirementTools::buildBackwardTree(backwardTree, this->mImposedScalarVariables, this->mImposedVectorVariables);

      stage.done();

      // Extract the run options for the equation parameters
      std::map<std::size_t,NonDimensional::SharedINumber> runOptions(this->mspEqParams->map());

      this->mspImposedTransformCoordinator = std::make_shared<Transform::TransformCoordinatorType>();

      // Initialise the transform coordinator
      std::vector<ArrayI> packs;
      Transform::TransformCoordinatorTools::computePacks(packs, this->mspImposedFwdGrouper, this->mspImposedBwdGrouper, {{0,forwardTree}}, {{0,backwardTree}}, {0}, this->mspRes);
      Transform::TransformCoordinatorTools::init(*this->mspImposedTransformCoordinator, this->mspImposedFwdGrouper, this->mspImposedBwdGrouper, packs, this->mspRes, runOptions);
   }

   void Coordinator::initSolvers()
   {
      StageTimer stage;
      stage.start("building trivial solvers");

      // Init trivial solver for trivial equations
      std::vector<SharedIScalarEquation> sEqs;
      std::vector<SharedIVectorEquation> vEqs;
      for(auto j: this->it())
      {
         auto sT = this->scalarRange(PseudospectralTag::Trivial::id(), j);
         sEqs.insert(sEqs.end(), sT.first, sT.second);
         auto vT = this->vectorRange(PseudospectralTag::Trivial::id(), j);
         vEqs.insert(vEqs.end(), vT.first, vT.second);
      }
      auto sT = std::make_pair(sEqs.begin(), sEqs.end());
      auto vT = std::make_pair(vEqs.begin(), vEqs.end());
      this->mTrivialCoordinator.init(sT, vT);

      stage.done();
      stage.start("building diagnostic solvers");

      // Init linear solver for trivial equations
      sEqs.clear();
      vEqs.clear();
      for(auto j: this->it())
      {
         auto sD = this->scalarRange(PseudospectralTag::Diagnostic::id(), j);
         sEqs.insert(sEqs.end(), sD.first, sD.second);
         auto vD = this->vectorRange(PseudospectralTag::Diagnostic::id(), j);
         vEqs.insert(vEqs.end(), vD.first, vD.second);
      }
      auto sD = std::make_pair(sEqs.begin(), sEqs.end());
      auto vD = std::make_pair(vEqs.begin(), vEqs.end());
      this->mLinearCoordinator.init(sD, vD);
      stage.done();
   }

   void Coordinator::updateEquationTime(const MHDFloat time, const bool finished)
   {
      for(auto j: this->it())
      {
         // Loop over all scalar equations
         assert(this->mScalarEquations.count(j) == 1);
         for(auto scalEqIt = this->mScalarEquations.at(j).begin(); scalEqIt < this->mScalarEquations.at(j).end(); ++scalEqIt)
         {
            (*scalEqIt)->setTime(time, finished);
         }

         // Loop over all vector equations
         assert(this->mVectorEquations.count(j) == 1);
         for(auto vectEqIt = this->mVectorEquations.at(j).begin(); vectEqIt < this->mVectorEquations.at(j).end(); ++vectEqIt)
         {
            (*vectEqIt)->setTime(time, finished);
         }
      }
   }

   void Coordinator::updateEquations(const int it, const bool isFinished) const
   {
      // Loop over all scalar equations
      assert(this->mScalarEquations.count(it) == 1);
      for(auto scalEqIt = this->mScalarEquations.at(it).begin(); scalEqIt < this->mScalarEquations.at(it).end(); ++scalEqIt)
      {
         (*scalEqIt)->updateConstraintKernel(this->time(), this->timestep(), isFinished);
      }

      // Loop over all vector equations
      assert(this->mVectorEquations.count(it) == 1);
      for(auto vectEqIt = this->mVectorEquations.at(it).begin(); vectEqIt < this->mVectorEquations.at(it).end(); ++vectEqIt)
      {
         (*vectEqIt)->updateConstraintKernel(this->time(), this->timestep(), isFinished);
      }
   }

   void Coordinator::writeDiagnostics(const bool isAsciiTime, const bool isHdf5Time) const
   {
      for(auto j: this->it())
      {
         // Loop over all scalar equations
         assert(this->mScalarEquations.count(j) == 1);
         for(auto scalEqIt = this->mScalarEquations.at(j).begin(); scalEqIt < this->mScalarEquations.at(j).end(); ++scalEqIt)
         {
            (*scalEqIt)->writeDiagnostics(isAsciiTime, isHdf5Time);
         }

         // Loop over all vector equations
         assert(this->mVectorEquations.count(j) == 1);
         for(auto vectEqIt = this->mVectorEquations.at(j).begin(); vectEqIt < this->mVectorEquations.at(j).end(); ++vectEqIt)
         {
            (*vectEqIt)->writeDiagnostics(isAsciiTime, isHdf5Time);
         }
      }
   }

   void Coordinator::computeNonlinear(const int it)
   {
      Profiler::RegionFixture fix("computeNonlinear");

      // Compute backward transform
      this->updatePhysical(it);

      // compute nonlinear interaction and forward transform
      this->updateSpectral(it);
   }

   void Coordinator::explicitTrivialEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      Profiler::RegionFixture<2> fix("Trivial-explicit");

      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Explicit term for trivial equations for operator " + ModelOperator::Coordinator::tag(opId), 5);
         this->mTrivialCoordinator.getExplicitInput(opId, scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      }
   }

   void Coordinator::explicitTrivialEquations(const std::size_t opId, const int it)
   {
      auto sT = this->scalarRange(PseudospectralTag::Trivial::id(), it);
      auto vT = this->vectorRange(PseudospectralTag::Trivial::id(), it);
      this->explicitTrivialEquations(opId, sT, vT);
   }

   void Coordinator::explicitDiagnosticEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      Profiler::RegionFixture<2> fix("Diagnostic-explicit");

      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Explicit term for diagnostic equations for operator " + ModelOperator::Coordinator::tag(opId), 5);
         this->mLinearCoordinator.getExplicitInput(opId, scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      }
   }

   void Coordinator::explicitDiagnosticEquations(const std::size_t opId, const int it)
   {
      auto sD = this->scalarRange(PseudospectralTag::Diagnostic::id(), it);
      auto vD = this->vectorRange(PseudospectralTag::Diagnostic::id(), it);
      this->explicitDiagnosticEquations(opId, sD, vD);
   }

   void Coordinator::explicitPrognosticEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      Profiler::RegionFixture<2> fix("Prognostic-explicit");

      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Explicit term for prognostic equations for operator " + ModelOperator::Coordinator::tag(opId), 5);
         this->mTimestepCoordinator.getExplicitInput(opId, scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      }
   }

   void Coordinator::explicitPrognosticEquations(const std::size_t opId, const int it)
   {
      auto sP = this->scalarRange(PseudospectralTag::Prognostic::id(), it);
      auto vP = this->vectorRange(PseudospectralTag::Prognostic::id(), it);
      this->explicitPrognosticEquations(opId, sP, vP);
   }

   void Coordinator::solveTrivialEquations(const std::size_t timeId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      Profiler::RegionFixture<2> fix("Trivial-solve");

      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Solve trivial equations at time " + SolveTiming::Coordinator::tag(timeId), 5);
         this->mTrivialCoordinator.setSolveTime(timeId);
         this->mTrivialCoordinator.solve(scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      }
   }

   void Coordinator::solveTrivialEquations(const std::size_t timeId, const int it)
   {
      auto sT = this->scalarRange(PseudospectralTag::Trivial::id(), it);
      auto vT = this->vectorRange(PseudospectralTag::Trivial::id(), it);
      this->solveTrivialEquations(timeId, sT, vT);
   }

   void Coordinator::solveDiagnosticEquations(const std::size_t timeId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      Profiler::RegionFixture<2> fix("Diagnostic-solve");

      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Solve diagnostic equations at time " + SolveTiming::Coordinator::tag(timeId), 5);
         this->mLinearCoordinator.setSolveTime(timeId);
         this->mLinearCoordinator.solve(scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      }
   }

   void Coordinator::solveDiagnosticEquations(const std::size_t timeId, const int it)
   {
      auto sD = this->scalarRange(PseudospectralTag::Diagnostic::id(), it);
      auto vD = this->vectorRange(PseudospectralTag::Diagnostic::id(), it);
      this->solveDiagnosticEquations(timeId, sD, vD);
   }

   void Coordinator::solvePrognosticEquations(ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      Profiler::RegionFixture fix("Prognostic-solve");

      // Only step forward if equations are present
      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Timestep prognostic equations", 5);
         this->mTimestepCoordinator.setSolveTime(SolveTiming::Prognostic::id());
         this->mTimestepCoordinator.stepForward(scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      }
   }

   void Coordinator::solvePrognosticEquations(const int it)
   {
      auto sP = this->scalarRange(PseudospectralTag::Prognostic::id(), it);
      auto vP = this->vectorRange(PseudospectralTag::Prognostic::id(), it);
      this->solvePrognosticEquations(sP, vP);
   }

   void Coordinator::explicitEquations(const int it)
   {
      Profiler::RegionFixture fix("explicitEquations");

      // Explicit trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitLinear::id(), it);

      // Explicit diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitLinear::id(), it);

      // Explicit prognostic equations
      this->explicitPrognosticEquations(ModelOperator::ExplicitLinear::id(), it);
   }

   void Coordinator::preSolveEquations()
   {
#warning "THIS NEEDS TO BE CHECK"
      StageTimer stage;
      stage.start("initializing fields");

      // only execute for for iteration
      // SHOULD THIS BE ITERATIONS UNTIL FIRST PROGNOSTIC EQUATIONS???
      int it = 0;

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveDiagnosticEquations(SolveTiming::After::id(), it);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveTrivialEquations(SolveTiming::After::id(), it);

      // Compute physical values
      this->updatePhysical(it);

      // Only compute forward transform for diagnostic and trivial equations
      this->updateSpectral(true, true, false, false, it);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solveDiagnosticEquations(SolveTiming::Before::id(), it);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solveTrivialEquations(SolveTiming::Before::id(), it);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveDiagnosticEquations(SolveTiming::After::id(), it);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveTrivialEquations(SolveTiming::After::id(), it);

      stage.done();

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

   void Coordinator::solveEquations(const int it)
   {
      Profiler::RegionFixture fix("solveEquations");

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solveTrivialEquations(SolveTiming::Before::id(), it);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solveDiagnosticEquations(SolveTiming::Before::id(), it);

      // Solve prognostic equations (timestep)
      this->explicitPrognosticEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solvePrognosticEquations(it);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveDiagnosticEquations(SolveTiming::After::id(), it);

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveTrivialEquations(SolveTiming::After::id(), it);
   }

   void Coordinator::finalizeTimestep()
   {
      // Update conditions at the end of timestep
      Profiler::RegionStart ("Control");
      if(this->mTimestepCoordinator.finishedStep())
      {
         // Update timestepper
         this->mTimestepCoordinator.update();

         // Update CFL condition
         this->mDiagnostics.updateCfl();

         // Synchronise diagnostics
         this->mDiagnostics.synchronize();

         // Adapt timestepper time step
         std::vector<SharedIScalarEquation> sEqs;
         std::vector<SharedIVectorEquation> vEqs;
         for(auto j: this->it())
         {
            auto sP = this->scalarRange(PseudospectralTag::Prognostic::id(), j);
            sEqs.insert(sEqs.end(), sP.first, sP.second);
            auto vP = this->vectorRange(PseudospectralTag::Prognostic::id(), j);
            vEqs.insert(vEqs.end(), vP.first, vP.second);
         }
         auto sP = std::make_pair(sEqs.begin(), sEqs.end());
         auto vP = std::make_pair(vEqs.begin(), vEqs.end());
         this->mTimestepCoordinator.adaptTimestep(this->mDiagnostics.cfl(), sP, vP);
      }
      Profiler::RegionStop ("Control");
   }

   void Coordinator::setupEquations()
   {
      // Get mesh from transform coordinator
      auto spMesh = std::make_shared<std::vector<Array> >(this->mTransformCoordinator.mesh());

      for(auto j: this->it())
      {
         // Loop over all scalar equations
         assert(this->mScalarEquations.count(j) == 1);
         for(auto scalEqIt = this->mScalarEquations.at(j).begin(); scalEqIt < this->mScalarEquations.at(j).end(); ++scalEqIt)
         {
            (*scalEqIt)->initSpectralMatrices();
            (*scalEqIt)->initConstraintKernel(spMesh);
            (*scalEqIt)->initSrcKernel();
         }

         // Loop over all vector equations
         assert(this->mVectorEquations.count(j) == 1);
         for(auto vectEqIt = this->mVectorEquations.at(j).begin(); vectEqIt < this->mVectorEquations.at(j).end(); ++vectEqIt)
         {
            (*vectEqIt)->initSpectralMatrices();
            (*vectEqIt)->initConstraintKernel(spMesh);
            (*vectEqIt)->initSrcKernel();
         }
      }

      // Share information across equations
      for(auto k: this->it())
      {
         for(auto sEqIt = this->mScalarEquations.at(k).begin(); sEqIt < this->mScalarEquations.at(k).end(); ++sEqIt)
         {
            for(auto j: this->it())
            {
               // Loop over all scalar equations
               assert(this->mScalarEquations.count(j) == 1);
               for(auto scalEqIt = this->mScalarEquations.at(j).begin(); scalEqIt < this->mScalarEquations.at(j).end(); ++scalEqIt)
               {
                  (*sEqIt)->linkEquation(*scalEqIt);
               }

               // Loop over all vector equations
               assert(this->mVectorEquations.count(j) == 1);
               for(auto vectEqIt = this->mVectorEquations.at(j).begin(); vectEqIt < this->mVectorEquations.at(j).end(); ++vectEqIt)
               {
                  (*sEqIt)->linkEquation(*vectEqIt);
               }
            }
         }

         for(auto vEqIt = this->mVectorEquations.at(k).begin(); vEqIt < this->mVectorEquations.at(k).end(); ++vEqIt)
         {
            for(auto j: this->it())
            {
               // Loop over all scalar equations
               assert(this->mScalarEquations.count(j) == 1);
               for(auto scalEqIt = this->mScalarEquations.at(j).begin(); scalEqIt < this->mScalarEquations.at(j).end(); ++scalEqIt)
               {
                  (*vEqIt)->linkEquation(*scalEqIt);
               }

               // Loop over all vector equations
               assert(this->mVectorEquations.count(j) == 1);
               for(auto vectEqIt = this->mVectorEquations.at(j).begin(); vectEqIt < this->mVectorEquations.at(j).end(); ++vectEqIt)
               {
                  (*vEqIt)->linkEquation(*vectEqIt);
               }
            }
         }
      }
   }

   void Coordinator::addToMap(const std::size_t eqId, const int it, const ScalarEquation_range& r)
   {
      auto key = std::make_pair(eqId, it);

      // Initialize map
      if(this->mScalarEqMap.count(key) == 0)
      {
         this->mScalarEqMap.insert(std::make_pair(key, std::vector<SharedIScalarEquation>()));
      }

      if(this->mScalarEqMap.at(key).size() > 0)
      {
         throw std::logic_error("Manually added Prognostic equation not supported");
      }
      else
      {
         auto ukey = std::make_pair(PseudospectralTag::Uninitialized::id(), it);

         for(auto eqIt = r.first; eqIt != r.second; eqIt++)
         {
            this->mScalarEqMap.at(key).push_back(*eqIt);
            auto ur = this->scalarRange(PseudospectralTag::Uninitialized::id(), it);
            auto uit = std::find(ur.first, ur.second, *eqIt);
            if(uit != ur.second)
            {
               this->mScalarEqMap.at(ukey).erase(uit);
            }
         }
      }
   }

   void Coordinator::addToMap(const std::size_t eqId, const int it, const VectorEquation_range& r)
   {
      auto key = std::make_pair(eqId, it);

      // Initialize map
      if(this->mVectorEqMap.count(key) == 0)
      {
         this->mVectorEqMap.insert(std::make_pair(key, std::vector<SharedIVectorEquation>()));
      }

      if(this->mVectorEqMap.at(key).size() > 0)
      {
         throw std::logic_error("Manually added Prognostic equation not supported");
      }
      else
      {
         auto ukey = std::make_pair(PseudospectralTag::Uninitialized::id(), it);

         for(auto eqIt = r.first; eqIt != r.second; eqIt++)
         {
            this->mVectorEqMap.at(key).push_back(*eqIt);
            auto ur = this->vectorRange(PseudospectralTag::Uninitialized::id(), it);
            auto uit = std::find(ur.first, ur.second, *eqIt);
            if(uit != ur.second)
            {
               this->mVectorEqMap.at(ukey).erase(uit);
            }
         }
      }
   }

   void Coordinator::sortEquations()
   {
      std::pair<int,int> pStart = {0,0};
      std::pair<int,int> dStart = {0,0};
      std::pair<int,int> tStart = {0,0};
      for(auto j: this->it())
      {
         // Sort scalar equations
         ScalarEquation_range sP;
         ScalarEquation_range sD;
         ScalarEquation_range sT;
         ScalarEquation_range sW;
         assert(this->mScalarEquations.count(j) == 1);
         Equations::Tools::sortByType(this->mScalarEquations.at(j), sP, sD, sT, sW);

         // Current implementation ony works if prognostic equations are in last iteration
         if(j != *this->it().rbegin() && std::distance(sP.first, sP.second) > 0)
         {
            throw std::logic_error("Current implementation requires scalar Prognostic equations to be in last sub-iteration");
         }

         // Add to map
         this->addToMap(PseudospectralTag::Prognostic::id(), j, sP);
         this->addToMap(PseudospectralTag::Diagnostic::id(), j, sD);
         this->addToMap(PseudospectralTag::Trivial::id(), j, sT);
         this->addToMap(PseudospectralTag::Wrapper::id(), j, sW);

         // Sort vector equations
         VectorEquation_range vP;
         VectorEquation_range vD;
         VectorEquation_range vT;
         VectorEquation_range vW;
         assert(this->mVectorEquations.count(j) == 1);
         Equations::Tools::sortByType(this->mVectorEquations.at(j), vP, vD, vT, vW);

         // Current implementation ony works if prognostic equations are in last iteration
         if(j != *this->it().rbegin() && std::distance(vP.first, vP.second) > 0)
         {
            throw std::logic_error("Current implementation requires vector Prognostic equations to be in last sub-iteration");
         }

         // Add to map
         this->addToMap(PseudospectralTag::Prognostic::id(), j, vP);
         this->addToMap(PseudospectralTag::Diagnostic::id(), j, vD);
         this->addToMap(PseudospectralTag::Trivial::id(), j, vT);
         this->addToMap(PseudospectralTag::Wrapper::id(), j, vW);

         // Identifiy the solver indexes by analysing the coupling between the equations
         DebuggerMacro_msg("Identifying solver for Prognostic equations sub-iteration " + std::to_string(j), 1);
         pStart = Equations::Tools::identifySolver(sP, vP, pStart.first, pStart.second);
         DebuggerMacro_msg("... done", 1);
         DebuggerMacro_msg("Identifying solver for Diagnostic equations sub-iteration " + std::to_string(j), 1);
         dStart = Equations::Tools::identifySolver(sD, vD, dStart.first, dStart.second);
         DebuggerMacro_msg("... done", 1);
         DebuggerMacro_msg("Identifying solver for Trivial equations sub-iteration " + std::to_string(j), 1);
         tStart = Equations::Tools::identifySolver(sT, vT, tStart.first, tStart.second);
         DebuggerMacro_msg("... done", 1);
      }
   }

   bool Coordinator::atLeastOne(ScalarEquation_range s, VectorEquation_range v) const
   {
      bool r = ((std::distance(s.first, s.second) + std::distance(v.first, v.second)) > 0);

      return r;
   }

   MHDFloat Coordinator::requiredStorage() const
   {
      MHDFloat mem = 0.0;
      #ifdef QUICC_STORAGEPROFILE
      #endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void Coordinator::profileStorage() const
   {
      #ifdef QUICC_STORAGEPROFILE
      // Profiling storage requirements
      this->mTransformCoordinator.profileStorage();
      #endif // QUICC_STORAGEPROFILE
   }
} // Pseudospectral
} // QuICC
