/**
 * @file Coordinator.cpp
 * @brief Source of the high level pseudospectral coordinator
 */

// System includes
//
#include <algorithm>
#include <stdexcept>
#include <type_traits>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Environment/QuICCEnv.hpp"
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
#include "QuICC/PseudospectralTag/Diagnostic.hpp"
#include "QuICC/PseudospectralTag/Prognostic.hpp"
#include "QuICC/PseudospectralTag/Trivial.hpp"
#include "QuICC/PseudospectralTag/Uninitialized.hpp"
#include "QuICC/PseudospectralTag/Wrapper.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "View/View.hpp"
#include "View/ViewUtils.hpp"
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


   namespace details {

   struct ptrAndIdxBlock
   {
      Memory::MemBlock<std::uint32_t> ptr;
      Memory::MemBlock<std::uint32_t> idx;
   };

   ptrAndIdxBlock getMeta(const TransformResolution& res, const std::uint32_t maxLayers, std::shared_ptr<Memory::memory_resource> mem)
   {
      std ::uint32_t nLayers = res.dim<QuICC::Dimensions::Data::DAT3D>();

      // ptr
      Memory::MemBlock<std::uint32_t> ptrBlock(maxLayers+1, mem.get());
      View::ViewBase<std::uint32_t> ptr(ptrBlock.data(), ptrBlock.size());

      std::uint32_t cumLayerSize = 0;
      std::uint32_t layerCounter = 0;
      ptr[0] = 0;
      for(std::uint32_t l = 0; l < maxLayers; ++l)
      {
         std::uint32_t layerSize = 0;
         if (layerCounter < nLayers)
         {
            auto layerIndex = static_cast<std::uint32_t>(res.idx<QuICC::Dimensions::Data::DAT3D>(layerCounter));
            if (l == layerIndex)
            {
               layerSize = res.dim<QuICC::Dimensions::Data::DAT2D>(layerCounter);
               ++layerCounter;
            }
         }
         ptr[l+1] = ptr[l] + layerSize;
         cumLayerSize += layerSize;
      }

      // idx
      Memory::MemBlock<std::uint32_t> idxBlock(cumLayerSize, mem.get());
      View::ViewBase<std::uint32_t> idx(idxBlock.data(), idxBlock.size());
      std::uint32_t l = 0;
      for(std::uint32_t k = 0; k < nLayers; ++k)
      {
         for(int j = 0; j < res.dim<QuICC::Dimensions::Data::DAT2D>(k); ++j)
         {
            // auto layerIndex = res.idx<QuICC::Dimensions::Data::DAT3D>(k);
            auto columnIndex = res.idx<QuICC::Dimensions::Data::DAT2D>(j, k);
            // auto columnHeight = res.dim<QuICC::Dimensions::Data::DATB1D>(j, k);
            idx[l] = columnIndex;
            ++l;
         }
      }

      ptrAndIdxBlock ret;
      ret.ptr = std::move(ptrBlock);
      ret.idx = std::move(idxBlock);
      return ret;
   }

   } // namespace details

   void Coordinator::addGraph(const std::string& graphStr)
   {
      /// Setup, get info from mspRes
      /// jw : look at BenchmarkMagC2.cpp:65
      /// LoadSplitter test?
      /// Transform::SharedTransformSetup WLFlBuilder::spSetup3D
      // mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // get Dims from mspRes
      // std::uint32_t Nr = jwRes.dim<Dimensions::Data::DATF1D>();
      std::uint32_t Nr = mspRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL);
      std::uint32_t N = mspRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // std::uint32_t Ntheta = alRes.dim<Dimensions::Data::DATF1D>();
      std::uint32_t Ntheta = mspRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::PHYSICAL);
      std::uint32_t L = mspRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      // std::uint32_t Nphi = ftRes.dim<Dimensions::Data::DATF1D>();
      std::uint32_t Nphi = mspRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::PHYSICAL);
      std::uint32_t M = mspRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Memory resource, depends on backend
      mMemRsr = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

      // get meta from mspRes
      const auto& jwRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA1D);
      auto metaJW = details::getMeta(jwRes, L, mMemRsr);
      const auto& alRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA2D);
      auto metaAL = details::getMeta(alRes, M, mMemRsr);
      const auto& ftRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA3D);
      auto metaFT = details::getMeta(ftRes, Nr, mMemRsr);

      constexpr std::uint32_t dim = 3;
      /// RThetaPhi - v012
      std::array<std::uint32_t, dim> physDims{Nr, Ntheta, Nphi};
      /// @brief Spectral dimensions
      /// NLM - v012
      std::array<std::uint32_t, dim> modsDims{N, L, M};

      // Layouts, depends on backend
      std::array<std::array<std::string, 2>, 3> layOpt;
      layOpt[0] = {"DCCSC3D", "DCCSC3D"};
      layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
      layOpt[2] = {"DCCSC3D", "DCCSC3D"};

      // Store meta stages to pass to Jitter
      std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
      meta.push_back({metaFT.ptr.data(), metaFT.ptr.size()});
      meta.push_back({metaFT.idx.data(), metaFT.idx.size()});
      meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
      meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
      meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
      meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

      mJitter = std::make_unique<QuICC::Graph::Jit<3>>(graphStr, mMemRsr, physDims, modsDims, layOpt, Graph::Stage::MMM, Graph::Stage::MMM, meta);

      // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
      std::array<View::ViewBase<std::uint32_t>, dim> pointersMods;
      pointersMods[1] = View::ViewBase<std::uint32_t>(metaJW.ptr.data(), metaJW.ptr.size());
      std::array<View::ViewBase<std::uint32_t>, dim> indicesMods;
      indicesMods[1] = View::ViewBase<std::uint32_t>(metaJW.idx.data(), metaJW.idx.size());

      // View for outputs/inputs
      std::vector<size_t> fields = {PhysicalNames::Temperature::id(), FieldComponents::Spectral::TOR, FieldComponents::Spectral::POL};

      // Add Views and Storage for each component
      auto scalarVarPtr = mspRes->sim().ss().bwdPtr(Dimensions::Transform::SPECTRAL);

      std::visit(
         [&](auto&& p)
         {
            // Get field scalar type
            using fld_t = typename std::remove_reference_t<decltype(p->data())>::Scalar;

            for (size_t f = 0; f < fields.size(); ++f)
            {
               // Field Id
               auto fId = fields[f];

               // Host mem block
               Memory::MemBlock<fld_t> block(modsDims[0]*metaJW.idx.size(), mMemRsr.get());

               // Host view
               // for now this works only for JW space
               std::array<std::uint32_t, dim> dims {modsDims[0], modsDims[2], modsDims[1]};
               View::View<fld_t, View::DCCSC3D> view(block.data(), block.size(), dims.data(), pointersMods.data(), indicesMods.data());

               // Store block
               mBlocksData.push_back(std::move(block));

               // Store view
               mId2View[fId] = view;
            }
         }, scalarVarPtr);

      // Store meta blocks
      mBlocksMeta.push_back(std::move(metaFT.ptr));
      mBlocksMeta.push_back(std::move(metaFT.idx));
      mBlocksMeta.push_back(std::move(metaAL.ptr));
      mBlocksMeta.push_back(std::move(metaAL.idx));
      mBlocksMeta.push_back(std::move(metaJW.ptr));
      mBlocksMeta.push_back(std::move(metaJW.idx));
   };

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

   void Coordinator::evolveBefore(const int curJ)
   {
      // Compute explicit linear terms
      DebuggerMacro_msg("Explicit equations", 4);
      this->explicitEquations(curJ);

      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::RUN);
      QuICCTimer().start();

      // Compute the nonlinear terms
      DebuggerMacro_msg("Transform loop", 4);
      this->computeNonlinear(curJ);

      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::NONLINEAR);
      QuICCTimer().update(ExecutionTimer::RUN);
      QuICCTimer().start();

      // Solve equations before prognostic
      this->solveEquationsBefore(curJ);
   }

   void Coordinator::evolveEndIteration(const int tIt)
   {
      // Solve equations after prognostic
      this->solveEquationsAfter(tIt);

      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::TIMESTEP);
      QuICCTimer().update(ExecutionTimer::RUN);
      QuICCTimer().start();
   }

   void Coordinator::evolveAfterPrognostic(const bool finishedStep)
   {
      auto tIt = *this->it().rbegin();

      this->evolveEndIteration(tIt);

      // Update the equations
      for(auto j: this->it())
      {
         DebuggerMacro_msg("Update equations", 4);
         this->updateEquations(j, finishedStep);
      }
   }

   void Coordinator::evolveUntilPrognostic(const bool finishedStep)
   {
      // Update equation time
      this->updateEquationTime(this->time(), finishedStep);

      // Loop over sub-steps
      for(auto j: this->it())
      {
         DebuggerMacro_showValue("Equation sub-iteration ", 3, j);

         // Evolve before prognostic
         this->evolveBefore(j);

         if(this->hasPrognostic(j))
         {
            // Solve prognostic equations (timestep)
            Profiler::RegionStart<2> ("Pseudospectral::Coordinator::solveEquations-prognostic");
            this->explicitPrognosticEquations(ModelOperator::ExplicitNonlinear::id(), j);
            Profiler::RegionStop<2> ("Pseudospectral::Coordinator::solveEquations-prognostic");
            break;
         }

         // End iteration
         this->evolveEndIteration(j);
      }
   }

   void Coordinator::evolve()
   {
      DebuggerMacro_msg("Evolving equations", 1);
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::evolve");

      // Solve prognostic equations (timestep)
      auto tIt = *this->it().rbegin();
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::solveEquations-prognostic");
      this->solvePrognosticEquations(tIt);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::solveEquations-prognostic");

      // Finish timestep iteration
      this->finalizeTimestep();
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
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::updatePhysical");

      // Define transform trees
      this->mTransformCoordinator.defineBwdTransforms(this->mBwdTree.at(it));

      // Compute physical values
      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mTransformCoordinator);
   }

   void Coordinator::updateSpectral(const int it)
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::updateSpectral");

      // Define transform trees
      this->mTransformCoordinator.defineFwdTransforms(this->mFwdTree.at(it));

      this->mspFwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mPhysicalKernels.at(it), this->mTransformCoordinator);
   }

   void Coordinator::updateSpectral(const bool isTrivial, const bool isDiagnostic, const bool isPrognostic, const bool isWrapper, const int it)
   {
      /// \todo This needs to be checked as it currently doesn't do anything

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
      this->mTimestepCoordinator.init(schemeId, this->mDiagnostics.startTime(), this->mDiagnostics.cfl(), this->mDiagnostics.maxError(), sP, vP, *this);

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
      const std::string profRegion = "Pseudospectral::Coordinator::updateEquations";
      Profiler::RegionFixture<1> fix(profRegion);

      // Loop over all scalar equations
      assert(this->mScalarEquations.count(it) == 1);
      for(auto scalEqIt = this->mScalarEquations.at(it).begin(); scalEqIt < this->mScalarEquations.at(it).end(); ++scalEqIt)
      {
         Profiler::RegionFixture<2> fixId(profRegion+"-"+PhysicalNames::Coordinator::tag((*scalEqIt)->name()));
         (*scalEqIt)->updateConstraintKernel(this->time(), this->timestep(), isFinished);
      }

      // Loop over all vector equations
      assert(this->mVectorEquations.count(it) == 1);
      for(auto vectEqIt = this->mVectorEquations.at(it).begin(); vectEqIt < this->mVectorEquations.at(it).end(); ++vectEqIt)
      {
         Profiler::RegionFixture<2> fixId(profRegion+"-"+PhysicalNames::Coordinator::tag((*vectEqIt)->name()));
         (*vectEqIt)->updateConstraintKernel(this->time(), this->timestep(), isFinished);
      }
   }

   void Coordinator::writeDiagnostics(const bool isAsciiTime, const bool isHdf5Time) const
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::writeDiagnostics");

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

namespace details
{
   template <class SCALAROUT, class SCALARIN>
   void copyEig2View(QuICC::View::View<SCALAROUT, View::DCCSC3D> view, const Eigen::Matrix<SCALARIN, -1, -1>& eig, const TransformResolution& res)
   {
      throw std::logic_error("trying to copy different types");
   }

   template <class SCALAR>
   void copyEig2View(QuICC::View::View<SCALAR, View::DCCSC3D> view, const Eigen::Matrix<SCALAR, -1, -1>& eig, const TransformResolution& res)
   {
      std ::uint32_t nLayers = res.dim<QuICC::Dimensions::Data::DAT3D>();

      /// copy data to view
      int start = 0;
      std::int64_t offSet = 0;
      for(std::uint32_t p = 0; p < nLayers; ++p)
      {
         // layer width
         int cols = res.dim<QuICC::Dimensions::Data::DAT2D>(p);
         // layer height
         int inRows = res.dim<QuICC::Dimensions::Data::DATB1D>(0, p);

         const Eigen::Ref<const Eigen::Matrix<SCALAR, -1, -1>> inB = eig.block(0, start, inRows, cols);

         for (std::int64_t j = 0; j < inB.cols(); ++j)
         {
            for (std::int64_t i = 0; i < inB.rows(); ++i)
            {
               #ifdef QUICC_JW_ROW_MAJOR
               // copy padded to flattened and transpose
               view[offSet + i*inB.cols()+j] = inB.data()[i+j*eig.rows()];
               #else
               // copy padded to flattened column
               view[offSet + i+j*inB.rows()] = inB.data()[i+j*eig.rows()];
               #endif
            }
         }
         offSet += inB.size();
         start += cols;
      }
   }

   template <class SCALAROUT, class SCALARIN>
   void copyView2Eig(Eigen::Matrix<SCALAROUT, -1, -1>& eig, const QuICC::View::View<SCALARIN, View::DCCSC3D> view, const TransformResolution& res)
   {
      throw std::logic_error("trying to copy different types");
   }

   template <class SCALAR>
   void copyView2Eig(Eigen::Matrix<SCALAR, -1, -1>& eig, const QuICC::View::View<SCALAR, View::DCCSC3D> view, const TransformResolution& res)
   {
      std ::uint32_t nLayers = res.dim<QuICC::Dimensions::Data::DAT3D>();

      int start = 0;
      std::uint64_t offSet = 0;
      for(std::uint32_t p = 0; p < nLayers; p++)
      {
         // layer width
         int cols = res.dim<QuICC::Dimensions::Data::DAT2D>(p);
         // layer height
         int outRows = res.dim<QuICC::Dimensions::Data::DATB1D>(0, p);

         Eigen::Ref<Eigen::Matrix<SCALAR, -1, -1>> outB = eig.block(0, start, outRows, cols);

         #ifdef QUICC_JW_ROW_MAJOR
         for (std::int64_t j = 0; j < outB.cols(); ++j)
         {
            for (std::int64_t i = 0; i < outB.rows(); ++i)
            {
               // copy padded to flattened column and transpose
               outB.data()[i+j*eig.rows()]= view[offSet + i*outB.cols()+j];
            }
         }
         #else
         for (std::int64_t i = 0; i < outB.size(); ++i)
         {
            outB.data()[i] = view[offSet + i];
         }
         #endif

         offSet += outB.size();
         start += cols;
      }
   }

} // namespace details

   void Coordinator::computeNonlinear(const int it)
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::computeNonlinear");

      assert(mJitter.get() != nullptr);

      // Copy to view
      const auto& jwRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA1D);
      auto temp = mScalarVariables[PhysicalNames::Temperature::id()];
      auto tempVarv = mId2View[PhysicalNames::Temperature::id()];
      std::visit(
         [&](auto&& p, auto& Tv)
         {
            auto ptrTemp = p->rDom(0).perturbation();
            details::copyEig2View(Tv, ptrTemp.data(), jwRes);
         }, temp, tempVarv);

      auto vec = mVectorVariables[PhysicalNames::Velocity::id()];
      auto TorVarv = mId2View[FieldComponents::Spectral::TOR];
      auto PolVarv = mId2View[FieldComponents::Spectral::POL];
      std::visit(
         [&](auto&& p, auto& Torv, auto& Polv)
         {
            auto ptrTor = p->rDom(0).perturbation().comp(FieldComponents::Spectral::TOR);
            details::copyEig2View(Torv, ptrTor.data(), jwRes);
            auto ptrPol = p->rDom(0).perturbation().comp(FieldComponents::Spectral::POL);
            details::copyEig2View(Polv, ptrPol.data(), jwRes);
         }, vec, TorVarv, PolVarv);


      // Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlOld");
      // // Compute backward transform
      // this->updatePhysical(it);

      // // compute nonlinear interaction and forward transform
      // this->updateSpectral(it);
      // Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlOld");

      Profiler::RegionStart<2>("Pseudospectral::Coordinator::nlNew");
      // Call graph
      std::visit(
         [&](auto& Tv, auto& Torv, auto& Polv)
         {
            mJitter->apply(Tv, Torv, Polv, Tv, Torv, Polv);
         }, tempVarv, TorVarv, PolVarv);
      Profiler::RegionStop<2>("Pseudospectral::Coordinator::nlNew");

      // Copy back
      std::visit(
         [&](auto&& p, auto& Tv)
         {
            auto ptrTemp = p->rDom(0).perturbation();
            details::copyView2Eig(ptrTemp.rData(), Tv, jwRes);
         }, temp, tempVarv);

      std::visit(
         [&](auto&& p, auto& Torv, auto& Polv)
         {
            auto ptrTor = p->rDom(0).perturbation().comp(FieldComponents::Spectral::TOR);
            details::copyView2Eig(ptrTor.rData(), Torv, jwRes);
            auto ptrPol = p->rDom(0).perturbation().comp(FieldComponents::Spectral::POL);
            details::copyView2Eig(ptrPol.rData(), Polv, jwRes);
         }, vec, TorVarv, PolVarv);
   }

   void Coordinator::explicitTrivialEquations(const std::size_t opId, ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
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

   bool Coordinator::hasPrognostic(const int it) const
   {
      auto sP = const_cast<Coordinator *>(this)->scalarRange(PseudospectralTag::Prognostic::id(), it);
      auto vP = const_cast<Coordinator *>(this)->vectorRange(PseudospectralTag::Prognostic::id(), it);
      return this->atLeastOne(sP, vP);
   }

   void Coordinator::solvePrognosticEquations(ScalarEquation_range scalarEq_range, VectorEquation_range vectorEq_range)
   {
      // Only step forward if equations are present
      if(this->atLeastOne(scalarEq_range, vectorEq_range))
      {
         DebuggerMacro_msg("Timestep prognostic equations", 5);
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
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::explicitEquations");

      // Explicit trivial equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::explicitEquations-trivial");
      this->explicitTrivialEquations(ModelOperator::ExplicitLinear::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::explicitEquations-trivial");

      // Explicit diagnostic equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::explicitEquations-diagnostic");
      this->explicitDiagnosticEquations(ModelOperator::ExplicitLinear::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::explicitEquations-diagnostic");

      // Explicit prognostic equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::explicitEquations-prognostic");
      this->explicitPrognosticEquations(ModelOperator::ExplicitLinear::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::explicitEquations-prognostic");
   }

   void Coordinator::preSolveEquations()
   {
      /// \todo No models are currently requiring a preSolve stage (only some cartesian models did).
      /// Implementation needs to be checked

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

   void Coordinator::solveEquationsBefore(const int it)
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::solveEquationsBefore");

      // Solve trivial equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::solveEquations-trivialBefore");
      this->explicitTrivialEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solveTrivialEquations(SolveTiming::Before::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::solveEquations-trivialBefore");

      // Solve diagnostic equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::solveEquations-diagnosticBefore");
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNonlinear::id(), it);
      this->solveDiagnosticEquations(SolveTiming::Before::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::solveEquations-diagnosticBefore");
   }

   void Coordinator::solveEquationsAfter(const int it)
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::solveEquationsAfter");

      // Solve diagnostic equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::solveEquations-diagnosticAfter");
      this->explicitDiagnosticEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveDiagnosticEquations(SolveTiming::After::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::solveEquations-diagnosticAfter");

      // Solve trivial equations
      Profiler::RegionStart<2> ("Pseudospectral::Coordinator::solveEquations-trivialAfter");
      this->explicitTrivialEquations(ModelOperator::ExplicitNextstep::id(), it);
      this->solveTrivialEquations(SolveTiming::After::id(), it);
      Profiler::RegionStop<2> ("Pseudospectral::Coordinator::solveEquations-trivialAfter");
   }

   void Coordinator::finalizeTimestep()
   {
      Profiler::RegionFixture<1> fix("Pseudospectral::Coordinator::finalizeTimestep");

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
