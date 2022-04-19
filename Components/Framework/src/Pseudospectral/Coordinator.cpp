/**
 * @file Coordinator.cpp
 * @brief Source of the high level pseudospectral coordinator
 */

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// Configuration includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <algorithm>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Pseudospectral/Coordinator.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
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
#include "QuICC/PseudospectralTag/Wrapper.hpp"

namespace QuICC {

namespace Pseudospectral {

   Coordinator::Coordinator()
   {
   }

   Coordinator::~Coordinator()
   {
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

   void Coordinator::addEquation(Equations::SharedIScalarEquation spEq)
   {
      // Set equation parameters if not yet set
      if(!this->mspEqParams)
      {
         this->mspEqParams = spEq->spEqParams();
      }

      this->mScalarEquations.push_back(spEq);
   }

   void Coordinator::addEquation(Equations::SharedIScalarEquation spEq, const std::size_t eqId, const bool addToList)
   {
      if(eqId == PseudospectralTag::Prognostic::id() ||
            eqId == PseudospectralTag::Diagnostic::id() ||
            eqId == PseudospectralTag::Trivial::id() ||
            eqId == PseudospectralTag::Wrapper::id())
      {
         throw std::logic_error("Prognostic, Diagnostic,Trivial and Wrapper equations are automatically added!");
      }

      if(addToList)
      {
         this->addEquation(spEq);
      }

      if(this->mScalarEqMap.count(eqId) == 0)
      {
         this->mScalarEqMap.insert(std::pair(eqId, std::vector<SharedIScalarEquation>()));
      }

      this->mScalarEqMap.find(eqId)->second.push_back(spEq);
   }

   void Coordinator::addEquation(Equations::SharedIVectorEquation spEq)
   {
      // Set equation parameters if not yet set
      if(!this->mspEqParams)
      {
         this->mspEqParams = spEq->spEqParams();
      }

      this->mVectorEquations.push_back(spEq);
   }

   void Coordinator::addEquation(Equations::SharedIVectorEquation spEq, const std::size_t eqId, const bool addToList)
   {
      if(eqId == PseudospectralTag::Prognostic::id() ||
            eqId == PseudospectralTag::Diagnostic::id() ||
            eqId == PseudospectralTag::Trivial::id() ||
            eqId == PseudospectralTag::Wrapper::id())
      {
         throw std::logic_error("Prognostic, Diagnostic,Trivial and Wrapper equations are automatically added!");
      }

      if(addToList)
      {
         this->addEquation(spEq);
      }

      if(this->mVectorEqMap.count(eqId) == 0)
      {
         this->mVectorEqMap.insert(std::pair(eqId, std::vector<SharedIVectorEquation>()));
      }

      this->mVectorEqMap.find(eqId)->second.push_back(spEq);
   }

   Coordinator::ScalarEquation_range Coordinator::scalarRange(const std::size_t eqId)
   {
      auto it = this->mScalarEqMap.find(eqId);
      auto r = std::make_pair(it->second.begin(), it->second.end());

      return r;
   }

   Coordinator::VectorEquation_range Coordinator::vectorRange(const std::size_t eqId)
   {
      auto it = this->mVectorEqMap.find(eqId);
      auto r = std::make_pair(it->second.begin(), it->second.end());

      return r;
   }

   Coordinator::SharedIScalarEquation Coordinator::scalarEq(const std::size_t eqId, const int i)
   {
      assert(this->mScalarEqMap.count(eqId) > 0);
      auto it = this->mScalarEqMap.find(eqId);

      assert(it->second.size() > 0);
      return it->second.at(i);
   }

   Coordinator::SharedIVectorEquation Coordinator::vectorEq(const std::size_t eqId, const int i)
   {
      assert(this->mVectorEqMap.count(eqId) > 0);
      auto it = this->mVectorEqMap.find(eqId);

      assert(it->second.size() > 0);
      return it->second.at(i);
   }

   void Coordinator::evolve()
   {
      bool isIntegrating = true;
      while(isIntegrating)
      {
         // Update equation time
         this->updateEquationTime(this->time(), this->mTimestepCoordinator.finishedStep());

         // Compute explicit linear terms
         this->explicitEquations();

         // Compute the nonlinear terms
         this->computeNonlinear();

         // Timestep the equations
         this->solveEquations();

         isIntegrating = !this->mTimestepCoordinator.finishedStep();
      }
   }

   void Coordinator::initTransformCoordinator(const std::vector<Transform::TransformTree>& integratorTree, const std::vector<Transform::TransformTree>& projectorTree)
   {
      // Extract the run options for the equation parameters
      std::map<std::size_t,NonDimensional::SharedINumber> runOptions(this->mspEqParams->map());

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(this->mTransformCoordinator, this->mspFwdGrouper, this->mspBwdGrouper, integratorTree, projectorTree, this->mspRes, runOptions);
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

   void Coordinator::updatePhysical()
   {
      // Compute physical values
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
   }

   void Coordinator::updateSpectral()
   {
      this->mspFwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mPhysicalKernels, this->mTransformCoordinator);
   }

   void Coordinator::updateSpectral(const bool isTrivial, const bool isDiagnostic, const bool isPrognostic, const bool isWrapper)
   {
      std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel> kernels;

      // Get kernels from trivial equations
      if(isTrivial)
      {
         auto sR = this->scalarRange(PseudospectralTag::Trivial::id());
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Trivial::id());
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Get kernels from diagnostic equations
      if(isDiagnostic)
      {
         auto sR = this->scalarRange(PseudospectralTag::Diagnostic::id());
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Diagnostic::id());
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Get kernels from prognostic equations
      if(isPrognostic)
      {
         auto sR = this->scalarRange(PseudospectralTag::Prognostic::id());
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Prognostic::id());
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Get kernels from wrapper equations
      if(isWrapper)
      {
         auto sR = this->scalarRange(PseudospectralTag::Wrapper::id());
         Equations::Tools::getNonlinearKernels(kernels, sR.first, sR.second);
         auto vR = this->vectorRange(PseudospectralTag::Wrapper::id());
         Equations::Tools::getNonlinearKernels(kernels, vR.first, vR.second);
      }

      // Set mesh for kernels
      Equations::Tools::setupPhysicalKernels(kernels, this->mTransformCoordinator.mesh());

      //this->mspFwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, kernels, this->mTransformCoordinator);
   }

   void Coordinator::init(const Array& tstep, const SharedSimulationBoundary spBcs)
   {
      StageTimer stage;

      stage.start("initializing variables");

      // Initialise the variables and set general variable requirements
      auto varInfo = RequirementTools::mergeRequirements(this->mScalarEquations, this->mVectorEquations);
      RequirementTools::initVariables(this->mScalarVariables, this->mVectorVariables, varInfo, this->mspRes);

      // Map variables to the equations and set nonlinear requirements
      std::vector<std::size_t> unmapped;
      RequirementTools::mapEquationVariables(this->mScalarEquations, this->mVectorEquations, this->mScalarVariables, this->mVectorVariables, spBcs, unmapped);

      if(unmapped.size() > 0)
      {
         throw std::logic_error("Variables not mapped to equations are currently not implemented");
      }

      // Initialize Imposed fields
      this->initImposed();

      // Transform trees
      std::vector<Transform::TransformTree> forwardTree;
      std::vector<Transform::TransformTree> backwardTree;
      RequirementTools::buildBackwardTree(backwardTree, this->mScalarEquations, this->mVectorEquations);
      RequirementTools::buildForwardTree(forwardTree, this->mScalarEquations, this->mVectorEquations);

      // Get nonlinear kernels
      Equations::Tools::getNonlinearKernels(this->mPhysicalKernels, this->mScalarEquations.begin(), this->mScalarEquations.end());
      Equations::Tools::getNonlinearKernels(this->mPhysicalKernels, this->mVectorEquations.begin(), this->mVectorEquations.end());

      stage.done();

      // Initialise the transform coordinator
      this->initTransformCoordinator(forwardTree, backwardTree);

      // Setup nonlinear kernels
      Equations::Tools::setupPhysicalKernels(this->mPhysicalKernels, this->mTransformCoordinator.mesh());

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

   void Coordinator::prepareEvolution()
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
      auto sP = this->scalarRange(PseudospectralTag::Prognostic::id());
      auto vP = this->vectorRange(PseudospectralTag::Prognostic::id());
      this->mTimestepCoordinator.init(this->mDiagnostics.startTime(), this->mDiagnostics.cfl(), this->mDiagnostics.maxError(), sP, vP);

      // Compute physical space values if required
      this->mspImposedBwdGrouper->transform(this->mImposedScalarVariables, this->mImposedVectorVariables, *this->mspImposedTransformCoordinator);
   }

   void Coordinator::cleanupForRun()
   {
      StageTimer stage;
      stage.start("Cleanup equation backends");

      // Loop over all scalar equations
      for(auto scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); ++scalEqIt)
      {
         (*scalEqIt)->cleanupBackend();
      }

      // Loop over all vector equations
      for(auto vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); ++vectEqIt)
      {
         (*vectEqIt)->cleanupBackend();
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
      auto varInfo = RequirementTools::mergeImposedRequirements(this->mScalarEquations, this->mVectorEquations);
      RequirementTools::initVariables(this->mImposedScalarVariables, this->mImposedVectorVariables, varInfo, this->mspRes);

      // Map variables to the equations and set nonlinear requirements
      std::vector<std::size_t> unmapped;
      RequirementTools::mapImposedVariables(this->mScalarEquations, this->mVectorEquations, this->mImposedScalarVariables, this->mImposedVectorVariables, unmapped);

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
      Transform::TransformCoordinatorTools::init(*this->mspImposedTransformCoordinator, this->mspImposedFwdGrouper, this->mspImposedBwdGrouper, forwardTree, backwardTree, this->mspRes, runOptions);
   }

   void Coordinator::initSolvers()
   {
      StageTimer stage;
      stage.start("building trivial solvers");

      // Init trivial solver for trivial equations
      auto sT = this->scalarRange(PseudospectralTag::Trivial::id());
      auto vT = this->vectorRange(PseudospectralTag::Trivial::id());
      this->mTrivialCoordinator.init(sT, vT);

      stage.done();
      stage.start("building diagnostic solvers");

      // Init linear solver for trivial equations
      auto sD = this->scalarRange(PseudospectralTag::Diagnostic::id());
      auto vD = this->vectorRange(PseudospectralTag::Diagnostic::id());
      this->mLinearCoordinator.init(sD, vD);
      stage.done();
   }

   void Coordinator::updateEquationTime(const MHDFloat time, const bool finished)
   {
      // Loop over all scalar equations
      for(auto scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); ++scalEqIt)
      {
         (*scalEqIt)->setTime(time, finished);
      }

      // Loop over all vector equations
      for(auto vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); ++vectEqIt)
      {
         (*vectEqIt)->setTime(time, finished);
      }
   }

   void Coordinator::computeNonlinear()
   {
      // Compute backward transform
      this->updatePhysical();

      // compute nonlinear interaction and forward transform
      this->mspFwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mPhysicalKernels, this->mTransformCoordinator);
   }

   void Coordinator::explicitTrivialEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      ProfilerMacro_start(Debug::Profiler::TRIVIALEQUATION);
      ProfilerMacro_start(Debug::Profiler::TRIVIALEXIN);
      this->mTrivialCoordinator.getExplicitInput(opId, scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(Debug::Profiler::TRIVIALEXIN);
      ProfilerMacro_stop(Debug::Profiler::TRIVIALEQUATION);
   }

   void Coordinator::explicitTrivialEquations(const std::size_t opId)
   {
      auto sT = this->scalarRange(PseudospectralTag::Trivial::id());
      auto vT = this->vectorRange(PseudospectralTag::Trivial::id());
      this->explicitTrivialEquations(opId, sT, vT);
   }

   void Coordinator::explicitDiagnosticEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      ProfilerMacro_start(Debug::Profiler::DIAGNOSTICEQUATION);
      ProfilerMacro_start(Debug::Profiler::DIAGNOSTICEXIN);
      this->mLinearCoordinator.getExplicitInput(opId, scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(Debug::Profiler::DIAGNOSTICEXIN);
      ProfilerMacro_stop(Debug::Profiler::DIAGNOSTICEQUATION);
   }

   void Coordinator::explicitDiagnosticEquations(const std::size_t opId)
   {
      auto sD = this->scalarRange(PseudospectralTag::Diagnostic::id());
      auto vD = this->vectorRange(PseudospectralTag::Diagnostic::id());
      this->explicitDiagnosticEquations(opId, sD, vD);
   }

   void Coordinator::explicitPrognosticEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      ProfilerMacro_start(Debug::Profiler::PROGNOSTICEQUATION);
      this->mTimestepCoordinator.getExplicitInput(opId, scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(Debug::Profiler::PROGNOSTICEQUATION);
   }

   void Coordinator::explicitPrognosticEquations(const std::size_t opId)
   {
      auto sP = this->scalarRange(PseudospectralTag::Prognostic::id());
      auto vP = this->vectorRange(PseudospectralTag::Prognostic::id());
      this->explicitPrognosticEquations(opId, sP, vP);
   }

   void Coordinator::solveTrivialEquations(const std::size_t timeId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      ProfilerMacro_start(Debug::Profiler::TRIVIALEQUATION);
      ProfilerMacro_start(Debug::Profiler::TRIVIALSOLVE);
      this->mTrivialCoordinator.setSolveTime(timeId);
      this->mTrivialCoordinator.solve(scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(Debug::Profiler::TRIVIALSOLVE);
      ProfilerMacro_stop(Debug::Profiler::TRIVIALEQUATION);
   }

   void Coordinator::solveTrivialEquations(const std::size_t timeId)
   {
      auto sT = this->scalarRange(PseudospectralTag::Trivial::id());
      auto vT = this->vectorRange(PseudospectralTag::Trivial::id());
      this->solveTrivialEquations(timeId, sT, vT);
   }

   void Coordinator::solveDiagnosticEquations(const std::size_t timeId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      ProfilerMacro_start(Debug::Profiler::DIAGNOSTICEQUATION);
      ProfilerMacro_start(Debug::Profiler::DIAGNOSTICSOLVE);
      this->mLinearCoordinator.setSolveTime(timeId);
      this->mLinearCoordinator.solve(scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(Debug::Profiler::DIAGNOSTICSOLVE);
      ProfilerMacro_stop(Debug::Profiler::DIAGNOSTICEQUATION);
   }

   void Coordinator::solveDiagnosticEquations(const std::size_t timeId)
   {
      auto sD = this->scalarRange(PseudospectralTag::Diagnostic::id());
      auto vD = this->vectorRange(PseudospectralTag::Diagnostic::id());
      this->solveDiagnosticEquations(timeId, sD, vD);
   }

   void Coordinator::solvePrognosticEquations(ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      ProfilerMacro_start(Debug::Profiler::PROGNOSTICEQUATION);
      this->mTimestepCoordinator.setSolveTime(SolveTiming::Prognostic::id());
      this->mTimestepCoordinator.stepForward(scalarEq_range, vectorEq_range, this->mScalarVariables, this->mVectorVariables);
      ProfilerMacro_stop(Debug::Profiler::PROGNOSTICEQUATION);
   }

   void Coordinator::solvePrognosticEquations()
   {
      auto sP = this->scalarRange(PseudospectralTag::Prognostic::id());
      auto vP = this->vectorRange(PseudospectralTag::Prognostic::id());
      this->solvePrognosticEquations(sP, vP);
   }

   void Coordinator::explicitEquations()
   {
      // Explicit trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitLinear::id());

      // Explicit diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitLinear::id());

      // Explicit prognostic equations
      this->explicitPrognosticEquations(ModelOperator::ExplicitLinear::id());
   }

   void Coordinator::preSolveEquations()
   {
      StageTimer stage;
      stage.start("initializing fields");

      /// \mhdBug This is not sufficient to recover all fields from previous computation

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id());
      this->solveDiagnosticEquations(SolveTiming::After::id());

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id());
      this->solveTrivialEquations(SolveTiming::After::id());

      // Compute physical values
      this->updatePhysical();

      // Only compute forward transform for diagnostic and trivial equations
      this->updateSpectral(true, true, false, false);

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNonlinear::id());
      this->solveDiagnosticEquations(SolveTiming::Before::id());

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNonlinear::id());
      this->solveTrivialEquations(SolveTiming::Before::id());

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id());
      this->solveDiagnosticEquations(SolveTiming::After::id());

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id());
      this->solveTrivialEquations(SolveTiming::After::id());

      stage.done();

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

   void Coordinator::solveEquations()
   {
      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNonlinear::id());
      this->solveTrivialEquations(SolveTiming::Before::id());

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNonlinear::id());
      this->solveDiagnosticEquations(SolveTiming::Before::id());

      // Solve prognostic equations (timestep)
      this->explicitPrognosticEquations(ModelOperator::ExplicitNonlinear::id());
      this->solvePrognosticEquations();

      // Solve diagnostic equations
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id());
      this->solveDiagnosticEquations(SolveTiming::After::id());

      // Solve trivial equations
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id());
      this->solveTrivialEquations(SolveTiming::After::id());

      // Update conditions at the end of timestep
      ProfilerMacro_start(Debug::Profiler::CONTROL);
      if(this->mTimestepCoordinator.finishedStep())
      {
         // Update timestepper
         this->mTimestepCoordinator.update();

         // Update CFL condition
         this->mDiagnostics.updateCfl();

         // Synchronise diagnostics
         this->mDiagnostics.synchronize();

         // Adapt timestepper time step
         auto sP = this->scalarRange(PseudospectralTag::Prognostic::id());
         auto vP = this->vectorRange(PseudospectralTag::Prognostic::id());
         this->mTimestepCoordinator.adaptTimestep(this->mDiagnostics.cfl(), sP, vP);
      }
      ProfilerMacro_stop(Debug::Profiler::CONTROL);
   }

   void Coordinator::setupEquations()
   {
      // Loop over all scalar equations
      for(auto scalEqIt = this->mScalarEquations.begin(); scalEqIt < this->mScalarEquations.end(); ++scalEqIt)
      {
         (*scalEqIt)->initSpectralMatrices();
         (*scalEqIt)->initConstraintKernel();
         (*scalEqIt)->initSrcKernel();
      }

      // Loop over all vector equations
      for(auto vectEqIt = this->mVectorEquations.begin(); vectEqIt < this->mVectorEquations.end(); ++vectEqIt)
      {
         (*vectEqIt)->initSpectralMatrices();
         (*vectEqIt)->initConstraintKernel();
         (*vectEqIt)->initSrcKernel();
      }
   }

   void Coordinator::sortEquations()
   {
      // Sort scalar equations
      ScalarEquation_range sP;
      ScalarEquation_range sD;
      ScalarEquation_range sT;
      ScalarEquation_range sW;
      Equations::Tools::sortByType(this->mScalarEquations, sP, sD, sT, sW);

      // Sort vector equations
      VectorEquation_range vP;
      VectorEquation_range vD;
      VectorEquation_range vT;
      VectorEquation_range vW;
      Equations::Tools::sortByType(this->mVectorEquations, vP, vD, vT, vW);

      // Identifiy the solver indexes by analysing the coupling between the equations
      Equations::Tools::identifySolver(sP, vP);
      Equations::Tools::identifySolver(sD, vD);
      Equations::Tools::identifySolver(sT, vT);

      // Add scalar equations to maps
      this->mScalarEqMap.insert(std::make_pair(PseudospectralTag::Prognostic::id(), std::vector<SharedIScalarEquation>()));
      auto sit = this->mScalarEqMap.find(PseudospectralTag::Prognostic::id());
      sit->second.insert(sit->second.end(), sP.first, sP.second);
      this->mScalarEqMap.insert(std::make_pair(PseudospectralTag::Diagnostic::id(), std::vector<SharedIScalarEquation>()));
      sit = this->mScalarEqMap.find(PseudospectralTag::Diagnostic::id());
      sit->second.insert(sit->second.end(), sD.first, sD.second);
      this->mScalarEqMap.insert(std::make_pair(PseudospectralTag::Trivial::id(), std::vector<SharedIScalarEquation>()));
      sit = this->mScalarEqMap.find(PseudospectralTag::Trivial::id());
      sit->second.insert(sit->second.end(), sT.first, sT.second);
      this->mScalarEqMap.insert(std::make_pair(PseudospectralTag::Wrapper::id(), std::vector<SharedIScalarEquation>()));
      sit = this->mScalarEqMap.find(PseudospectralTag::Wrapper::id());
      sit->second.insert(sit->second.end(), sW.first, sW.second);

      // Add vector equations to maps
      this->mVectorEqMap.insert(std::make_pair(PseudospectralTag::Prognostic::id(), std::vector<SharedIVectorEquation>()));
      auto vit = this->mVectorEqMap.find(PseudospectralTag::Prognostic::id());
      vit->second.insert(vit->second.end(), vP.first, vP.second);
      this->mVectorEqMap.insert(std::make_pair(PseudospectralTag::Diagnostic::id(), std::vector<SharedIVectorEquation>()));
      vit = this->mVectorEqMap.find(PseudospectralTag::Diagnostic::id());
      vit->second.insert(vit->second.end(), vD.first, vD.second);
      this->mVectorEqMap.insert(std::make_pair(PseudospectralTag::Trivial::id(), std::vector<SharedIVectorEquation>()));
      vit = this->mVectorEqMap.find(PseudospectralTag::Trivial::id());
      vit->second.insert(vit->second.end(), vT.first, vT.second);
      this->mVectorEqMap.insert(std::make_pair(PseudospectralTag::Wrapper::id(), std::vector<SharedIVectorEquation>()));
      vit = this->mVectorEqMap.find(PseudospectralTag::Wrapper::id());
      vit->second.insert(vit->second.end(), vW.first, vW.second);
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
}
}
