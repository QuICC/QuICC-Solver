/**
 * @file Interface.cpp
 * @brief Implementation of the timestep coordinator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Timestep/Interface.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Timestep/Constants.hpp"

namespace QuICC {

namespace Timestep {

   Interface::Interface(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
      : mcMinCnst(2), mcMaxJump(MAX_STEPSIZE_JUMP), mcUpWindow(FIVE_PC_WINDOW), mcMinDt(LIMIT_MINSTEP), mcMaxDt(LIMIT_MAXSTEP), mMaxError(-1.0), mOldDt(this->mcMinDt), mDt(2,1), mTime(0.0), mRefTime(0.0), mCnstSteps(0.0), mStepTime(0.0)
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

}
}
