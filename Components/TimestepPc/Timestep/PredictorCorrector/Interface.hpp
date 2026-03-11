/**
 * @file Interface.hpp
 * @brief Implementation of an interface to the ImEx predictor-corrector schemes
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACE_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Timestep/Constants.hpp"
#include "QuICC/Timestep/IScheme.hpp"
#include "QuICC/Timestep/Interface.hpp"
#include "QuICC/Timestep/TimestepSolverCoordinator.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

/**
 * @brief Implementation of general timestepper structure
 */
template <typename TScheme> class Interface : public Timestep::Interface
{
public:
   /// Typedef for solver implementation
   template <typename T1, typename T2, template <typename> class T3>
   using SolverImplementationType =
      typename TScheme::template ImplementationType<T1, T2, T3>;

   /// Typedef for parent coordinator
   typedef TimestepSolverCoordinator<SolverImplementationType>
      SolverCoordinator;

   /// Typedef for a shared real operator solver
   typedef
      typename SolverCoordinator::SharedRealSolverType SharedRealSolverType;

   /// Typedef for a shared complex operator solver
   typedef typename SolverCoordinator::SharedComplexSolverType
      SharedComplexSolverType;

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
    *
    * \mhdBug Not fully implemented
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
    * @brief Timestep is finished?
    */
   bool finishedStep() const;

   /**
    * @brief Set solve time
    */
   void setSolveTime(const std::size_t timeId);
private:
   /**
    * @brief Update time dependence
    */
   void updateMatrices();

   /**
    * @brief Interface to timestepping scheme
    */
   SharedIScheme mspScheme;

   /**
    * @brief Interface to timestepping scheme
    */
   SolverCoordinator mSolverCoord;
};

template <typename TScheme>
Interface<TScheme>::Interface(const MHDFloat time, const Matrix& cfl,
   const MHDFloat maxError, const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq, Pseudospectral::Coordinator& pseudo) :
    Timestep::Interface(time, cfl, maxError, scalEq, vectEq, pseudo)
{
   // Create Timestepper scheme
   std::shared_ptr<TScheme> spScheme = std::make_shared<TScheme>();

   // Use embedded scheme to compute error
   if (maxError > 0.0)
   {
      spScheme->enableEmbedded();
      this->mMaxError = maxError;
   }
   this->mspScheme = spScheme;

   this->mSolverCoord.init(this->timestep(), scalEq, vectEq, spScheme);
}

template <typename TScheme>
void Interface<TScheme>::tuneAdaptive(const MHDFloat time)
{
   this->mStepTime = time;
}

template <typename TScheme> bool Interface<TScheme>::finishedStep() const
{
   return this->mSolverCoord.finishedStep();
}

template <typename TScheme>
void Interface<TScheme>::setSolveTime(const std::size_t timeId)
{
   this->mSolverCoord.setSolveTime(timeId);
}

template <typename TScheme>
void Interface<TScheme>::getExplicitInput(const std::size_t opId,
   const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq,
   const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
{
   Profiler::RegionFixture<2> fix("Timestep-explicitInput");

   this->mSolverCoord.getExplicitInput(opId, scalEq, vectEq, scalVar, vectVar);
}

template <typename TScheme>
void Interface<TScheme>::adaptTimestep(const Matrix& cfl)
{
   // Process CFL information
   this->processCfl(cfl, this->mSolverCoord.error(), this->mspScheme->order());

   //
   // Update the timestep matrices if necessary
   //
   if (this->timestep() != this->mOldDt && this->timestep() > 0.0)
   {
      DebuggerMacro_showValue(
         "Updating timestep and matrices with new Dt = ", 0, this->timestep());

      this->mSolverCoord.updateTimestep(this->timestep());

      DebuggerMacro_start("Update matrices", 0);
      // Update the time dependence in matrices
      this->updateMatrices();
      DebuggerMacro_stop("Update matrices t = ", 0);

      DebuggerMacro_start("Complex operator update", 0);
      // Update solvers from complex operator, complex field steppers
      Solver::updateSolvers<SolverImplementationType,
         typename Solver::SparseCoordinatorBase<
            SolverImplementationType>::ComplexSolver_iterator>(
         this->mSolverCoord);
      DebuggerMacro_stop("Complex operator solver update t = ", 0);

      DebuggerMacro_start("Real operator solver update", 0);
      // Update solvers from real operator, complex field steppers
      Solver::updateSolvers<SolverImplementationType,
         typename Solver::SparseCoordinatorBase<
            SolverImplementationType>::RealSolver_iterator>(this->mSolverCoord);
      DebuggerMacro_stop("Real operator solver update t = ", 0);
   }
   else
   {
      this->mCnstSteps += 1.0;
   }

   // Update CFL writer
   this->writeCfl();
}

template <typename TScheme>
void Interface<TScheme>::stepForward(const ScalarEquation_range& scalEq,
   const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar,
   const VectorVariable_map& vectVar)
{
   bool isIntegrating = true;
   while (isIntegrating)
   {
      DebuggerMacro_msg("Time integration sub-step", 2);

      this->mpPseudo->evolveUntilPrognostic(this->finishedStep());
      this->setSolveTime(SolveTiming::Prognostic::id());

      Profiler::RegionStart<2>("Timestep-input");
      // Update the equation input to the timestepper
      this->mSolverCoord.getInput(scalEq, vectEq, scalVar, vectVar);
      Profiler::RegionStop<2>("Timestep-input");

      Profiler::RegionStart<2>("Timestep-solve");
      // Solve all the linear systems
      this->mSolverCoord.solveSystems();
      Profiler::RegionStop<2>("Timestep-solve");

      Profiler::RegionStart<2>("Timestep-output");
      // Transfer timestep output back to equations
      this->mSolverCoord.transferOutput(scalEq, vectEq);
      Profiler::RegionStop<2>("Timestep-output");

      // Clear the solver RHS
      this->mSolverCoord.clearSolvers();

      // Update current time
      this->mTime =
         this->mRefTime + this->mSolverCoord.stepFraction() * this->timestep();

      this->mpPseudo->evolveAfterPrognostic(this->finishedStep());

      isIntegrating = !this->finishedStep();
   }
}

template <typename TScheme> void Interface<TScheme>::updateMatrices()
{
   // Loop over all complex operator, complex field timesteppers
   Solver::updateTimeMatrixSolvers<SolverImplementationType,
      typename Solver::SparseCoordinatorBase<
         SolverImplementationType>::ComplexSolver_iterator>(this->mSolverCoord,
      this->timestep());

   // Loop over all real operator, complex field timesteppers
   Solver::updateTimeMatrixSolvers<SolverImplementationType,
      typename Solver::SparseCoordinatorBase<
         SolverImplementationType>::RealSolver_iterator>(this->mSolverCoord,
      this->timestep());
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

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_INTERFACE_HPP
