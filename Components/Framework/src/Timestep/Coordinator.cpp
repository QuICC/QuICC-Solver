/**
 * @file Coordinator.cpp
 * @brief Implementation of the timestep coordinator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Timestep/Coordinator.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Timestep {

   Coordinator::Coordinator()
      : ParentCoordinator(), mcMinCnst(2), mcMaxJump(1.602), mcUpWindow(1.05), mcMinDt(1e-10), mcMaxDt(1e-1), mMaxError(-1.0), mOldDt(this->mcMinDt), mDt(2,1), mTime(0.0), mRefTime(0.0), mCnstSteps(0.0), mStepTime(0.0)
   {
      this->mDt(0,0) = this->mcMinDt;
      this->mDt(1,0) = -100.0;

      // Create CFL writer
      auto spCflWriter = std::make_shared<Io::Ascii::CflWriter>();
      this->mspIo = spCflWriter;
      this->mspIo->init();
   }

   Coordinator::~Coordinator()
   {
      this->mspIo->finalize();
   }

   void Coordinator::tuneAdaptive(const MHDFloat time)
   {
      this->mStepTime = time;
   }

   MHDFloat Coordinator::time() const
   {
      return this->mTime;
   }

   MHDFloat Coordinator::timestep() const
   {
      return this->mDt(0,0);
   }

   void Coordinator::update()
   {
      this->mTime = this->mRefTime + this->timestep();
      this->mRefTime = this->mTime;
   }

   void Coordinator::adaptTimestep(const Matrix& cfl, const ScalarEquation_range&, const VectorEquation_range&)
   {
      // Store old timestep
      this->mOldDt = this->timestep();

      // Update CFL information
      this->mDt.block(0, 1, this->mDt.rows(), cfl.cols()-1) = cfl.rightCols(cfl.cols()-1);

      // New computed CFL
      MHDFloat compCfl = cfl(0,0);

      // Check if CFL allows for a larger timestep
      MHDFloat newCflDt = 0.0;
      if(compCfl > this->mcUpWindow*this->timestep())
      {
         if(this->mCnstSteps >= this->mcMinCnst)
         {
            // Set new timestep
            newCflDt = std::min(compCfl, this->mcMaxJump*this->timestep());
         } else
         {
            // Reuse same timestep
            newCflDt = this->timestep();
         }

      // Check if CFL is below minimal timestep or downard jump is large
      } else if(compCfl < this->mcMinDt || compCfl < this->timestep()/this->mcMaxJump)
      {
         // Signal simulation abort
         newCflDt = -compCfl;

      // Check if CFL requires a lower timestep
      } else if(compCfl < this->timestep()*(2.0-this->mcUpWindow))
      {
         // Set new timestep
         newCflDt = compCfl;

      } else
      {
         newCflDt = this->timestep();
      }

      // Gather error across processes
      #ifdef QUICC_MPI
      if(this->mError > 0.0)
      {
         MPI_Allreduce(MPI_IN_PLACE, &this->mError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      }
      #endif //QUICC_MPI

      // No error control and no CFL condition
      MHDFloat newErrorDt = 0.0;

      // Use what ever condition is used by CFL
      if(this->mError < 0)
      {
         newErrorDt = -1.0;

      // Error is too large, reduce timestep
      } else if(this->mError > this->mMaxError)
      {
         newErrorDt = this->timestep()*std::pow(this->mMaxError/this->mError,1./this->mspScheme->order())/this->mcUpWindow;

      // Error is small, increase timestep
      } else if(this->mError < this->mMaxError/(this->mcMaxJump*0.9) && this->mCnstSteps >= this->mcMinCnst)
      {
         newErrorDt = std::min(this->timestep()*std::pow(this->mMaxError/this->mError,1./this->mspScheme->order()), this->timestep()*this->mcMaxJump);

      // Timestep should not be increased
      } else
      {
         newErrorDt = this->timestep();
      }

      // Update error details
      if(this->mMaxError > 0.0)
      {
         this->mDt(0, this->mDt.cols()-1) = newErrorDt;
         this->mDt(1, this->mDt.cols()-1) = this->mError;
      }

      // CFL condition requested abort!
      if(newCflDt < 0.0)
      {
         this->mDt(0,0) = newCflDt;
         this->mDt(1,0) = cfl(1,0);

      // Get minimum between both conditions
      } else if(newCflDt > 0.0 && newErrorDt > 0.0)
      {
         if(newCflDt < newErrorDt)
         {
            this->mDt(0,0) = newCflDt;
            this->mDt(1,0) = cfl(1,0);
         } else
         {
            this->mDt(0,0) = newErrorDt;
            this->mDt(1,0) = -200.0;
         }

      // Use CFL condition
      } else if(newCflDt > 0.0)
      {
         if(this->timestep() != newCflDt)
         {
            this->mDt(0,0) = newCflDt;
            this->mDt(1,0) = cfl(1,0);
         }

      // Use error condition
      } else if(newErrorDt > 0.0)
      {
         this->mDt(0,0) = newErrorDt;
         this->mDt(1,0) = -200.0;
      }

      //
      // Update the timestep matrices if necessary
      //
      if(this->timestep() != this->mOldDt && this->timestep() > 0.0)
      {
         DebuggerMacro_showValue("Updating timestep and matrices with new Dt = ", 0, this->timestep());

         DebuggerMacro_start("Update matrices", 0);
         // Update the time dependence in matrices
         this->updateMatrices();
         DebuggerMacro_stop("Update matrices t = ", 0);

         DebuggerMacro_start("Complex operator update", 0);
         // Update solvers from complex operator, complex field steppers
         Solver::updateSolvers<TimeSchemeSelector::ImplementationType,typename Solver::SparseCoordinatorBase<TimeSchemeSelector::ImplementationType>::ComplexSolver_iterator>(*this);
         DebuggerMacro_stop("Complex operator solver update t = ", 0);

         DebuggerMacro_start("Real operator solver update", 0);
         // Update solvers from real operator, complex field steppers
         Solver::updateSolvers<TimeSchemeSelector::ImplementationType,typename Solver::SparseCoordinatorBase<TimeSchemeSelector::ImplementationType>::RealSolver_iterator>(*this);
         DebuggerMacro_stop("Real operator solver update t = ", 0);
      } else
      {
         this->mCnstSteps += 1.0;
      }

      // Update CFL writer
      this->mspIo->setSimTime(this->mTime, this->mDt, this->mCnstSteps);
      this->mspIo->write();

      if(this->timestep() != this->mOldDt && this->timestep() > 0.0)
      {
         this->mCnstSteps = 0.0;
      }
   }

   void Coordinator::stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      ProfilerMacro_start(Debug::Profiler::TSTEPIN);
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq, scalVar, vectVar);
      ProfilerMacro_stop(Debug::Profiler::TSTEPIN);

      ProfilerMacro_start(Debug::Profiler::TSTEPSOLVE);
      // Solve all the linear systems
      this->solveSystems();
      ProfilerMacro_stop(Debug::Profiler::TSTEPSOLVE);

      ProfilerMacro_start(Debug::Profiler::TSTEPOUT);
      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);
      ProfilerMacro_stop(Debug::Profiler::TSTEPOUT);

      // Clear the solver RHS
      this->clearSolvers();

      // Update current time
      this->mTime = this->mRefTime + this->stepFraction()*this->timestep();
   }

   void Coordinator::init(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      // Set initial time
      this->mTime = time;
      this->mRefTime = this->mTime;

      // Create Timestepper scheme
      std::shared_ptr<TimeSchemeSelector> spScheme = std::make_shared<TimeSchemeSelector>();
      spScheme->init();

      // Use embedded scheme to compute error
      if(maxError > 0.0)
      {
         spScheme->enableEmbedded();
         this->mMaxError = maxError;
      }
      this->mspScheme = spScheme;

      // Set initial timestep
      this->mOldDt = cfl(0,0);

      // Update CFL details
      if(this->mMaxError > 0.0)
      {
         this->mDt.resize(cfl.rows(), cfl.cols()+1);
         this->mDt.leftCols(cfl.cols()) = cfl;
         this->mDt.rightCols(1)(0) = 0.0;
         this->mDt.rightCols(1)(1) = -200.0;
      } else
      {
         this->mDt = cfl;
      }

      DebuggerMacro_showValue("Creating timestepper with initial timestep Dt = ", 0, this->timestep());

      // Create solvers
      ParentCoordinator::createCoordSolvers(scalEq, vectEq);

      // Set Timestepper scheme
      for(auto solIt = this->mRealSolvers.begin(); solIt != this->mRealSolvers.end(); ++ solIt)
      {
         (*solIt)->setScheme(spScheme);
      }

      for(auto solIt = this->mComplexSolvers.begin(); solIt != this->mComplexSolvers.end(); ++ solIt)
      {
         (*solIt)->setScheme(spScheme);
      }

      // Initialise solver storage
      ParentCoordinator::initCoordSolvers(scalEq, vectEq);
   }

   void Coordinator::updateMatrices()
   {
      // Loop over all complex operator, complex field timesteppers
      Solver::updateTimeMatrixSolvers<TimeSchemeSelector::ImplementationType,typename Solver::SparseCoordinatorBase<TimeSchemeSelector::ImplementationType>::ComplexSolver_iterator>(*this, this->timestep());

      // Loop over all real operator, complex field timesteppers
      Solver::updateTimeMatrixSolvers<TimeSchemeSelector::ImplementationType,typename Solver::SparseCoordinatorBase<TimeSchemeSelector::ImplementationType>::RealSolver_iterator>(*this, this->timestep());
   }

   void Coordinator::buildSolverMatrix(Coordinator::SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildTimestepMatrixWrapper(spSolver, spEq, comp, this->timestep(), idx);
   }

   void Coordinator::buildSolverMatrix(Coordinator::SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildTimestepMatrixWrapper(spSolver, spEq, comp, this->timestep(), idx);
   }

   void Coordinator::printInfo(std::ostream& stream)
   {
      // Create nice looking ouput header
      Tools::Formatter::printNewline(stream);
      Tools::Formatter::printLine(stream, '-');
      Tools::Formatter::printCentered(stream, "Timestepper information", '*');
      Tools::Formatter::printLine(stream, '-');

      std::stringstream oss;
      int base = 20;

      // Timestep scheme
      oss << "Timestepper: " << this->mspScheme->name() << " (" << this->mspScheme->order() << ")";
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
      #endif //defined QUICC_SPLINALG_MUMPS

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
      #endif //defined QUICC_SPTRILINALG_SPARSELU

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
      #endif //defined QUICC_SPSPDLINALG_SIMPLICIALLDLT

      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      Tools::Formatter::printLine(stream, '*');
      Tools::Formatter::printNewline(stream);
   }


}
}
