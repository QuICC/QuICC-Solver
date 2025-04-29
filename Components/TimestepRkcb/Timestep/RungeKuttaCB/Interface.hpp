/**
 * @file Interface.hpp
 * @brief Implementation of a interface to the RK CB schemes
 */

#ifndef QUICC_TIMESTEP_RUNGEKUTTACB_INTERFACE_HPP
#define QUICC_TIMESTEP_RUNGEKUTTACB_INTERFACE_HPP

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

namespace RungeKuttaCB {

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
    * @brief Timestep is finished?
    */
   bool finishedStep() const final;

   /**
    * @brief Set solve time
    */
   void setSolveTime(const std::size_t timeId) final;

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
    * @param scalEq  Shared scalar equations
    * @param vectEq  Shared vector equations
    */
   void adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq,
      const VectorEquation_range& vectEq) final;

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
   this->mSolverCoord.getExplicitInput(opId, scalEq, vectEq, scalVar, vectVar);
}

template <typename TScheme>
void Interface<TScheme>::adaptTimestep(const Matrix& cfl,
   const ScalarEquation_range&, const VectorEquation_range&)
{
   // Store old timestep
   this->mOldDt = this->timestep();

   // Update CFL information
   this->mDt.block(0, 1, this->mDt.rows(), cfl.cols() - 1) =
      cfl.rightCols(cfl.cols() - 1);

   // New computed CFL
   MHDFloat compCfl = cfl(0, 0);

   // Check if CFL allows for a larger timestep
   MHDFloat newCflDt = 0.0;
   if (compCfl > this->mcUpWindow * this->timestep())
   {
      if (this->mCnstSteps >= this->mcMinCnst)
      {
         // Set new timestep
         newCflDt = std::min(compCfl, this->mcMaxJump * this->timestep());
      }
      else
      {
         // Reuse same timestep
         newCflDt = this->timestep();
      }

      // Check if CFL is below minimal timestep or downard jump is large
   }
   else if (compCfl < this->mcMinDt ||
            compCfl < this->timestep() / this->mcMaxJump)
   {
      // Signal simulation abort
      newCflDt = -compCfl;

      // Check if CFL requires a lower timestep
   }
   else if (compCfl < this->timestep() * (2.0 - this->mcUpWindow))
   {
      // Set new timestep
      newCflDt = compCfl;
   }
   else
   {
      newCflDt = this->timestep();
   }

   // Get timestepper error (if applicable)
   MHDFloat error = this->mSolverCoord.error();

// Gather error across processes
#ifdef QUICC_MPI
   if (error > 0.0)
   {
      MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX,
         MPI_COMM_WORLD);
   }
#endif // QUICC_MPI

   // No error control and no CFL condition
   MHDFloat newErrorDt = 0.0;

   // Use what ever condition is used by CFL
   if (error < 0)
   {
      newErrorDt = -1.0;

      // Error is too large, reduce timestep
   }
   else if (error > this->mMaxError)
   {
      newErrorDt =
         this->timestep() *
         std::pow(this->mMaxError / error, 1. / this->mspScheme->order()) /
         this->mcUpWindow;

      // Error is small, increase timestep
   }
   else if (error < this->mMaxError / (this->mcMaxJump * 0.9) &&
            this->mCnstSteps >= this->mcMinCnst)
   {
      newErrorDt =
         std::min(this->timestep() * std::pow(this->mMaxError / error,
                                        1. / this->mspScheme->order()),
            this->timestep() * this->mcMaxJump);

      // Timestep should not be increased
   }
   else
   {
      newErrorDt = this->timestep();
   }

   // Update error details
   if (this->mMaxError > 0.0)
   {
      this->mDt(0, this->mDt.cols() - 1) = newErrorDt;
      this->mDt(1, this->mDt.cols() - 1) = error;
   }

   // CFL condition requested abort!
   if (newCflDt < 0.0)
   {
      this->mDt(0, 0) = newCflDt;
      this->mDt(1, 0) = cfl(1, 0);

      // Get minimum between both conditions
   }
   else if (newCflDt > 0.0 && newErrorDt > 0.0)
   {
      if (newCflDt < newErrorDt)
      {
         this->mDt(0, 0) = newCflDt;
         this->mDt(1, 0) = cfl(1, 0);
      }
      else
      {
         this->mDt(0, 0) = newErrorDt;
         this->mDt(1, 0) = ERROR_LOCATION;
      }

      // Use CFL condition
   }
   else if (newCflDt > 0.0)
   {
      if (this->timestep() != newCflDt)
      {
         this->mDt(0, 0) = newCflDt;
         this->mDt(1, 0) = cfl(1, 0);
      }

      // Use error condition
   }
   else if (newErrorDt > 0.0)
   {
      this->mDt(0, 0) = newErrorDt;
      this->mDt(1, 0) = ERROR_LOCATION;
   }

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
   this->mspIo->setSimTime(this->mTime, this->mDt, this->mCnstSteps);
   this->mspIo->write();

   if (this->timestep() != this->mOldDt && this->timestep() > 0.0)
   {
      this->mCnstSteps = 0.0;
   }
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
   // Create nice looking ouput header
   Tools::Formatter::printNewline(stream);
   Tools::Formatter::printLine(stream, '-');
   Tools::Formatter::printCentered(stream, "Timestepper information", '*');
   Tools::Formatter::printLine(stream, '-');

   std::stringstream oss;
   int base = 20;

   // Timestep scheme
   oss << "Timestepper: " << this->mspScheme->name() << " ("
       << this->mspScheme->order() << ")";
   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   // General linear solver
   oss << "General solver: ";
#if defined QUICC_SPLINALG_MUMPS
   oss << "MUMPS";
#elif defined QUICC_SPLINALG_UMFPACK
   oss << "UmfPack";
#elif defined QUICC_SPLINALG_SPARSELU
   oss << "SparseLU";
#else
   oss << "(unknown)";
#endif // defined QUICC_SPLINALG_MUMPS

   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   // Triangular linear solver
   oss << "Triangular solver: ";
#if defined QUICC_SPTRILINALG_SPARSELU
   oss << "SparseLU";
#elif defined QUICC_SPTRILINALG_MUMPS
   oss << "MUMPS";
#elif defined QUICC_SPTRILINALG_UMFPACK
   oss << "UmfPack";
#else
   oss << "(unknown)";
#endif // defined QUICC_SPTRILINALG_SPARSELU

   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   // SPD linear solver
   oss << "SPD solver: ";
#if defined QUICC_SPSPDLINALG_SIMPLICIALLDLT
   oss << "SimplicialLDLT";
#elif defined QUICC_SPSPDLINALG_SIMPLICIALLLT
   oss << "SimplicialLLT";
#elif defined QUICC_SPSPDLINALG_MUMPS
   oss << "MUMPS";
#elif defined QUICC_SPSPDLINALG_UMFPACK
   oss << "UmfPack";
#elif defined QUICC_SPSPDLINALG_SPARSELU
   oss << "SparseLU";
#else
   oss << "(unknown)";
#endif // defined QUICC_SPSPDLINALG_SIMPLICIALLDLT

   Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
   oss.str("");

   Tools::Formatter::printLine(stream, '*');
   Tools::Formatter::printNewline(stream);
}

} // namespace RungeKuttaCB
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_RUNGEKUTTACB_INTERFACE_HPP
