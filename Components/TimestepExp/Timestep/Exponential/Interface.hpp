/**
 * @file Interface.hpp
 * @brief Implementation of an interface to the exponential schemes
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_INTERFACE_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_INTERFACE_HPP

// System includes
//
#include <memory>
#include <type_traits>

// Project includes
//
#include "Memory/MemoryResource.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/PhysicalNames/JacobianMagnetic.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/PseudospectralTag/Prognostic.hpp"
#include "QuICC/Register/Rhs.hpp"
#include "QuICC/Register/Solution.hpp"
#include "QuICC/Register/Temporary.hpp"
#include "QuICC/Timestep/IScheme.hpp"
#include "QuICC/Timestep/Interface.hpp"
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "Timestep/Exponential/TimestepperCoordinator.hpp"
#include "Timestep/Exponential/EpirkTimestepper.hpp"
#include "Timestep/Exponential/InterfaceFunctors.hpp"
#include "Timestep/Exponential/AugmentedJacobianFunctor.hpp"
#include "Timestep/Exponential/Functors/CallExplicitPrognosticFunctor.hpp"
#include "Memory/Cpu/NewDelete.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Implementation of general interface structure for exponential schemes
 */
template <typename TScheme> class Interface : public Timestep::Interface
{
public:
   /// Typedef for solver implementation
   template <typename T1, typename T2, typename T3>
   using SolverImplementationType =
      EpirkTimestepper<T1, T2, T3>;

   /// Typedef for parent coordinator
   typedef TimestepperCoordinator<SolverImplementationType, base_t>
      SolverCoordinator;
   /**
    * @brief Constructor
    *
    * @param time    Initial time value
    * @param cfl     Initial CFL timestep
    * @param error   Max error allowed during timestep
    * @param scalEq  Shared scalar equations
    * @param vectEq  Shared vector equations
    */
   Interface(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
      const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
      Pseudospectral::Coordinator& pseudo);

   /**
    * @brief Destructor
    */
   virtual ~Interface() = default;

   /**
    * @brief Update equation explicit linear input to solver
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void getExplicitInput(const std::size_t opId,
      const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
      const ScalarVariable_map& scalVar,
      const VectorVariable_map& vectVar) final;

   /**
    * @brief Tune adaptive timestepper
    */
   void tuneAdaptive(const MHDFloat time) final;

   /**
    * @brief Adapt the timestep used
    *
    * @param cfl     CFL conditions
    */
   void adaptTimestep(const Matrix& cfl) final;

   /**
    * @brief Compute (partial) forward step
    *
    * @param scalEq Shared scalar equations
    * @param vectEq Shared vector equations
    * @param scalVar Shared scalar variables
    * @param vectVar Shared vector variables
    */
   void stepForward(const ScalarEquation_range& scalEq,
      const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar,
      const VectorVariable_map& vectVar) final;

   /**
    * @brief Print timestepper information to stream
    *
    * @param stream  Output stream
    */
   void printInfo(std::ostream& stream) final;

protected:
   using Timestep::Interface::printInfo;

   /**
    * @brief Initialize functors
    */
   void initFunctors();

   /**
    * @brief Initialize solution
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

   /**
    * @brief Update equation input to solver
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

   /**
    * @brief Transfer solution from solver
    *
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

   /**
    * @brief Translate equations to timestepper info
    *
    * @param infos  Vector of information
    * @param scalEq Scalar equations
    * @param vectEq Vector equations
    */
   void translate(std::vector<TimestepperInfo>& infos,
      const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

private:
   std::set<int>  mBaseItIds;

   /**
    * @brief Interface to timestepping scheme
    */
   SharedIScheme mspScheme;

   /**
    * @brief Interface to timestepper coordinator
    */
   SolverCoordinator mSolverCoord;

   /**
    * @brief
    */
   std::shared_ptr<Memory::memory_resource> _mem;

   /**
    * @brief Shared field ID to solver field id
    */
   std::shared_ptr<std::map<SpectralFieldId, std::size_t>> mpFieldIdMap;

   /**
    * @brief Shared Jacobian functor
    */
   std::shared_ptr<AugmentedJacobianFunctor> mpJac;

   using TSFunctor = GetStepperFunctor<SolverCoordinator>;
   std::shared_ptr<TSFunctor> mpTsFunc;
};

template <typename TScheme>
Interface<TScheme>::Interface(const MHDFloat time, const Matrix& cfl,
   const MHDFloat maxError, const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq, Pseudospectral::Coordinator& pseudo) :
    Timestep::Interface(time, cfl, maxError, scalEq, vectEq, pseudo), mBaseItIds({0})
{
   _mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();

   // Create Timestepper scheme
   std::shared_ptr<TScheme> spScheme = std::make_shared<TScheme>();

   // Use embedded scheme to compute error
   if (maxError > 0.0)
   {
      this->mMaxError = maxError;
   }
   this->mspScheme = spScheme;

   // Init functors
   this->initFunctors();

   std::vector<TimestepperInfo> infos;
   this->translate(infos, scalEq, vectEq);

   assert(this->mpFieldIdMap);
   this->mpJac = std::make_shared<AugmentedJacobianFunctor>(-1.0, Register::Temporary::id(), 0, 1, this->mpPseudo, this->mpFieldIdMap, _mem);

   this->mSolverCoord.init(this->timestep(), infos, spScheme, this->mpJac);

   this->initSolution(scalEq, vectEq);
}

template <typename TScheme>
void Interface<TScheme>::initFunctors()
{
   this->mpTsFunc = std::make_shared<TSFunctor>(this->mSolverCoord);
}

template <typename TScheme>
void Interface<TScheme>::translate(std::vector<TimestepperInfo>& infos, const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq)
{
   this->mpFieldIdMap = std::make_shared<std::map<SpectralFieldId, std::size_t>>();
   auto bFunc = std::make_shared<DoNothingFunctor>();
   auto pFunc = std::make_shared<TranslateInfoFunctor>(infos, this->mpFieldIdMap, this->_mem);
   auto aFunc = std::make_shared<DoNothingFunctor>();
   ProcessRangeFunctor processor(bFunc, pFunc, aFunc, *this->mBaseItIds.begin());
   processor(scalEq);
   processor(vectEq);

   // Update system size
   std::size_t sysN = 0;
   for(auto&& t: infos)
   {
      sysN += t.rows;
   }
   for(auto&& t: infos)
   {
      t.rows = sysN;
   }

   // Add Jacobians to field ID map
   std::map<SpectralFieldId, SpectralFieldId> jacMap = {
      {std::make_pair(PhysicalNames::JacobianTemperature::id(), FieldComponents::Spectral::SCALAR), 
         std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR)},
      {std::make_pair(PhysicalNames::JacobianVelocity::id(), FieldComponents::Spectral::TOR), 
         std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR)},
      {std::make_pair(PhysicalNames::JacobianVelocity::id(), FieldComponents::Spectral::POL), 
         std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL)},
      {std::make_pair(PhysicalNames::JacobianMagnetic::id(), FieldComponents::Spectral::TOR), 
         std::make_pair(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::TOR)},
      {std::make_pair(PhysicalNames::JacobianMagnetic::id(), FieldComponents::Spectral::POL), 
         std::make_pair(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::POL)}
   };
   for(auto&& [jId, id]: jacMap)
   {
      if(this->mpFieldIdMap->count(id) > 0)
      {
         this->mpFieldIdMap->emplace(jId, this->mpFieldIdMap->at(id));
      }
   }
}

template <typename TScheme>
void Interface<TScheme>::tuneAdaptive(const MHDFloat time)
{
   this->mStepTime = time;
}

template <typename TScheme>
void Interface<TScheme>::initSolution(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   DebuggerMacro_msg("Initialize timestepper solutions", 6);

   auto bFunc = std::make_shared<DoNothingFunctor>();
   using OpFunctor = InitSolutionFunctor<TSFunctor>;
   auto pvFunc = std::make_shared<OpFunctor>(this->mpTsFunc);
   auto pFunc = std::make_shared<InputFunctor<OpFunctor>>(pvFunc, this->mpFieldIdMap, this->_mem);
   auto aFunc = std::make_shared<DoNothingFunctor>();
   ProcessRangeFunctor processor(bFunc, pFunc, aFunc, *this->mBaseItIds.begin());
   processor(scalEq);
   processor(vectEq);
}

template <typename TScheme>
void Interface<TScheme>::getExplicitInput(const std::size_t opId,
   const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
   const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
{
   Profiler::RegionFixture<2> fix("Timestep-explicitInput");
   throw std::logic_error("This should not be called");
}

template <typename TScheme>
void Interface<TScheme>::getInput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   Profiler::RegionFixture<2> fix("Timestep-input");

   auto bFunc = std::make_shared<ApplyConstraintFunctor>(SolveTiming::Before::id());
   using OpFunctor = GetInputFunctor<TSFunctor>;
   auto pvFunc = std::make_shared<OpFunctor>(this->mpTsFunc, Register::Rhs::id(), 1);
   auto pFunc = std::make_shared<InputFunctor<OpFunctor>>(pvFunc, this->mpFieldIdMap, this->_mem);
   auto aFunc = std::make_shared<DoNothingFunctor>();
   ProcessRangeFunctor processor(bFunc, pFunc, aFunc, *this->mBaseItIds.begin());
   processor(scalEq);
   processor(vectEq);
}

template <typename TScheme>
void Interface<TScheme>::transferOutput(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
{
   Profiler::RegionFixture<2> fix("Timestep-output");

   auto bFunc = std::make_shared<DoNothingFunctor>();
   using OpFunctor = TransferOutputFunctor<TSFunctor>;
   auto pvFunc = std::make_shared<OpFunctor>(this->mpTsFunc, Register::Solution::id(), 0);
   using CorrFunctor = TransferCorrectionFunctor<TSFunctor>;
   auto cvFunc = std::make_shared<CorrFunctor>(this->mpTsFunc, Register::Solution::id(), 0);
   auto pFunc = std::make_shared<OutputFunctor<OpFunctor, CorrFunctor>>(pvFunc, cvFunc, this->mpFieldIdMap, this->_mem);
   auto aFunc = std::make_shared<DoNothingFunctor>();
   ProcessRangeFunctor processor(bFunc, pFunc, aFunc, *this->mBaseItIds.begin());
   processor(scalEq);
   processor(vectEq);
}

template <typename TScheme>
void Interface<TScheme>::adaptTimestep(const Matrix& cfl)
{
   // Process CFL information
   this->processCfl(cfl, -1, this->mspScheme->order());

   //
   // Update the timestep matrices if necessary
   //
   if (this->timestep() != this->mOldDt && this->timestep() > 0.0)
   {
      DebuggerMacro_showValue(
         "Updating timestep to new Dt = ", 0, this->timestep());

      // Update the time dependencies
      DebuggerMacro_start("Update time dependencies", 0);
      this->mSolverCoord.updateTimestep(this->timestep());
      DebuggerMacro_stop("Update time dependencies", 0);
   }
   else
   {
      this->mCnstSteps += 1.0;
   }

   // Update CFL writer
   this->writeCfl();
}

template <typename TScheme>
void Interface<TScheme>::stepForward(const ScalarEquation_range& scalEq_ignore,
   const VectorEquation_range& vectEq_ignore, const ScalarVariable_map& scalVar,
   const VectorVariable_map& vectVar)
{
   DebuggerMacro_msg("Time integration sub-step", 2);

   int curIt = *this->mBaseItIds.begin();
   auto scalEq = this->mpPseudo->scalarRange(PseudospectralTag::Prognostic::id(), curIt);
   auto vectEq = this->mpPseudo->vectorRange(PseudospectralTag::Prognostic::id(), curIt);

   auto progFunc = std::make_shared<Functors::CallExplicitPrognosticFunctor<TSFunctor>>(this->mpTsFunc, Register::Rhs::id(), 1, *this->mBaseItIds.begin(), this->mpFieldIdMap, this->_mem);
   this->mpPseudo->evolveUntilPrognostic(this->mBaseItIds, false, progFunc);

   // Update the equation input to the timestepper
   this->getInput(scalEq, vectEq);

   Profiler::RegionStart<2>("Timestep-solve");
   int jacIt = 1;
   auto jscalEq = this->mpPseudo->scalarRange(PseudospectralTag::Prognostic::id(), jacIt);
   auto jvectEq = this->mpPseudo->vectorRange(PseudospectralTag::Prognostic::id(), jacIt);
   this->mpJac->setEquations(jscalEq, jvectEq);
   // Solve all the linear systems
   this->mSolverCoord.stepForward();
   Profiler::RegionStop<2>("Timestep-solve");

   // Transfer timestep output back to equations
   this->transferOutput(scalEq, vectEq);

   // Reset solvers
   this->mSolverCoord.resetSolvers();

   // Update current time
   this->mTime =
      this->mRefTime + this->timestep();

   this->mpPseudo->evolveAfterPrognostic(this->mBaseItIds, true);
}

template <typename TScheme>
void Interface<TScheme>::printInfo(std::ostream& stream)
{
   // Timestep scheme
   std::stringstream oss;
   oss << "Timestepper: " << this->mspScheme->name() << " ("
       << this->mspScheme->order() << ")";

   this->printInfo(stream, oss.str());
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_INTERFACE_HPP
