/**
 * @file Interface.cpp
 * @brief Implementation of the timestep coordinator
 */

// System includes
//

// Project includes
//
#include "QuICC/Timestep/Interface.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Timestep/Constants.hpp"

namespace QuICC {

namespace Timestep {

   Interface::Interface(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, Pseudospectral::Coordinator& pseudo)
      : mcMinCnst(2), mcMaxJump(MAX_STEPSIZE_JUMP), mcUpWindow(FIVE_PC_WINDOW), mcMinDt(LIMIT_MINSTEP), mcMaxDt(LIMIT_MAXSTEP), mMaxError(-1.0), mOldDt(this->mcMinDt), mDt(2,1), mTime(0.0), mRefTime(0.0), mCnstSteps(0.0), mStepTime(0.0), mpPseudo(&pseudo)
   {
      this->mDt(0,0) = this->mcMinDt;
      this->mDt(1,0) = FIXEDSTEP_LOCATION;

      // Create CFL writer
      auto spCflWriter = std::make_shared<Io::Ascii::CflWriter>();
      this->mspIo = spCflWriter;
      this->mspIo->init();

      // Set initial time
      this->mTime = time;
      this->mRefTime = this->mTime;

      // Set initial timestep
      this->mOldDt = cfl(0,0);

      // Update CFL details
      if(this->mMaxError > 0.0)
      {
         this->mDt.resize(cfl.rows(), cfl.cols()+1);
         this->mDt.leftCols(cfl.cols()) = cfl;
         this->mDt.rightCols(1)(0) = 0.0;
         this->mDt.rightCols(1)(1) = ERROR_LOCATION;
      } else
      {
         this->mDt = cfl;
      }

      DebuggerMacro_showValue("Creating timestepper with initial timestep Dt = ", 0, this->timestep());
   }

   Interface::~Interface()
   {
      this->mspIo->finalize();
   }

   MHDFloat Interface::time() const
   {
      return this->mTime;
   }

   MHDFloat Interface::timestep() const
   {
      return this->mDt(0,0);
   }

   void Interface::update()
   {
      this->mTime = this->mRefTime + this->timestep();
      this->mRefTime = this->mTime;
   }

   void Interface::processCfl(const Matrix& cfl, const MHDFloat err, const int order)
   {
      // Copy to simplify MPI code
      MHDFloat error = err;

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
            std::pow(this->mMaxError / error, 1. / order) /
            this->mcUpWindow;

         // Error is small, increase timestep
      }
      else if (error < this->mMaxError / (this->mcMaxJump * 0.9) &&
            this->mCnstSteps >= this->mcMinCnst)
      {
         newErrorDt =
            std::min(this->timestep() * std::pow(this->mMaxError / error,
                     1. / order),
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
   }

   void Interface::writeCfl()
   {
      // Update CFL writer
      this->mspIo->setSimTime(this->mTime, this->mDt, this->mCnstSteps);
      this->mspIo->write();

      if (this->timestep() != this->mOldDt && this->timestep() > 0.0)
      {
         this->mCnstSteps = 0.0;
      }
   }

void Interface::printInfo(std::ostream& stream, const std::string schemeInfo)
{
   // Create nice looking ouput header
   Tools::Formatter::printNewline(stream);
   Tools::Formatter::printLine(stream, '-');
   Tools::Formatter::printCentered(stream, "Timestepper information", '*');
   Tools::Formatter::printLine(stream, '-');

   std::stringstream oss;
   int base = 20;

   // Timestep scheme
   oss << "Timestepper: " << schemeInfo;
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

}
}
