/** 
 * @file BreakPoint.cpp
 * @brief Source of the profiler breakpoint
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/BreakPoint.hpp"

// Project includes
//

#include <iostream>
namespace QuICC {

namespace Debug {

namespace Profiler {

   BreakPoint shiftPoint(BreakPoint point, const int shift)
   {
      return static_cast<BreakPoint>(static_cast<int>(point) + shift);
   }

   std::string pointName(BreakPoint point)
   {
      unsigned int val = static_cast<int>(point);
  
      // LVL 3 names for Worland transform
      if(static_cast<unsigned int>(WORLANDTRA) < val && val - static_cast<unsigned int>(WORLANDTRA) < LVL0)
      {
         switch(val % LVL2)
         {
            case 0:
               break;
            case 1:
               return "Pre";
            case 2:
               return "FFT";
            case 3:
               return "Post";
            default:
               return "ID = " + std::to_string(val % LVL2);
         }
      }

      switch(point)
      {
         case BLACKHOLE:
            return "Blackhole";

         //////////////////////////
         case BWDTRANSFORM:
            return "Backward";

         case BWDIN:
            return "Input";
         case BWD1D:
            return "1D";
         case BWD1DIN:
            return "In";
         case BWD1DTRA:
            return "Transform";
         case BWD1DTRA_A:
            return "A";
         case BWD1DTRA_B:
            return "B";
         case BWD1DTRA_C:
            return "C";
         case BWD1DTRA_D:
            return "D";
         case BWD1DTRA_E:
            return "E";
         case BWD1DTRA_F:
            return "F";
         case BWD1DTRA_G:
            return "G";
         case BWD1DTRA_H:
            return "H";
         case BWD1DTRA_I:
            return "I";
         case BWD1DOUT:
            return "Out";
         case BWD1DOUTWORK:
            return "Work";
         case BWD1DOUTCOMM:
            return "Comm";
         case BWD2D:
            return "2D";
         case BWD2DIN:
            return "In";
         case BWD2DTRA:
            return "Transform";
         case BWD2DOUT:
            return "Out";
         case BWD2DOUTWORK:
            return "Work";
         case BWD2DOUTCOMM:
            return "Comm";
         case BWDND:
            return "ND";
         case BWDNDIN:
            return "In";
         case BWDNDTRA:
            return "Transform";
         case BWDNDOUT:
            return "Out";
         case BWDNDOUTWORK:
            return "Work";
         case BWDNDOUTCOMM:
            return "Comm";
         case BWDOUT:
            return "Output";

         //////////////////////////
         case NONLINEAR:
            return "Nonlinear";

         //////////////////////////
         case FWDTRANSFORM: 
            return "Forward";

         case FWDIN:
            return "Input";
         case FWDND:
            return "ND";
         case FWDNDIN:
            return "In";
         case FWDNDTRA:
            return "Transform";
         case FWDNDOUT:
            return "Out";
         case FWDNDOUTWORK:
            return "Work";
         case FWDNDOUTCOMM:
            return "Comm";
         case FWD2D:
            return "2D";
         case FWD2DIN:
            return "In";
         case FWD2DTRA:
            return "Transform";
         case FWD2DOUT:
            return "Out";
         case FWD2DOUTWORK:
            return "Work";
         case FWD2DOUTCOMM:
            return "Comm";
         case FWD1D:
            return "1D";
         case FWD1DIN:
            return "In";
         case FWD1DTRA:
            return "Transform";
         case FWD1DOUT:
            return "Out";
         case FWD1DOUTWORK:
            return "Work";
         case FWD1DOUTCOMM:
            return "Comm";
         case FWDOUT:
            return "Output";

         //////////////////////////
         case PROGNOSTICEQUATION:
            return "Prognostic";

         case TSTEPEXIN:
            return "Input(Ex)";
         case TSTEPIN:
            return "Input";
         case TSTEPSOLVE:
            return "Solve";
         case TSTEPOUT:
            return "Output";

         //////////////////////////
         case DIAGNOSTICEQUATION:
            return "Diagnostic";

         case DIAGNOSTICSOLVE:
            return "Solve";
         case DIAGNOSTICEXIN:
            return "Input(Ex)";

         //////////////////////////
         case TRIVIALEQUATION:
            return "Trivial";

         case TRIVIALSOLVE:
            return "Solve";
         case TRIVIALEXIN:
            return "Input(Ex)";

         //////////////////////////
         case CONTROL:
            return "Control";

         //////////////////////////
         case IO:
            return "IO";

         //////////////////////////
         case CONSTRAIN:
            return "Constrain";

         //////////////////////////
         case WORLANDTRA:
            return "Worland transform";

         case WORLANDPROJ:
            return "Projector";
         case WORLANDPROJ_P:
            return "P";
         case WORLANDPROJ_D1:
            return "D1";
         case WORLANDPROJ_D1R1:
            return "D1R1";
         case WORLANDPROJ_DIVR1D1R1:
            return "DivR1D1";
         case WORLANDPROJ_DIVR1D1R1_ZERO:
            return "DivR1D1_Zero";
         case WORLANDPROJ_DIVR1:
            return "DivR1";
         case WORLANDPROJ_DIVR1_ZERO:
            return "DivR1_Zero";
         case WORLANDPROJ_SPHLAPL:
            return "SphLapl";

         case WORLANDINTG:
            return "Integrator";
         case WORLANDINTG_P:
            return "P";
         case WORLANDINTG_R1:
            return "R1";
         case WORLANDINTG_DIVR1:
            return "DivR1";
         case WORLANDINTG_DIVR1D1R1:
            return "DivR1D1R1";
         case WORLANDINTG_I2:
            return "I2";
         case WORLANDINTG_I2DIVR1D1R1_ZERO:
            return "I2DivR1D1R1_Zero";
         case WORLANDINTG_I2DIVR1_ZERO:
            return "I2DivR1_Zero";
         case WORLANDINTG_I4DIVR1D1R1_ZERO:
            return "I4DivR1D1R1_Zero";
         case WORLANDINTG_I4DIVR1_ZERO:
            return "I4DivR1_Zero";
         case WORLANDINTG_P_ZERO:
            return "P_Zero";
         case WORLANDINTG_I2_ZERO:
            return "I2_Zero";
         case WORLANDINTG_DIVR1_ZERO:
            return "DivR1_Zero";
         case WORLANDINTG_DIVR1D1R1_ZERO:
            return "DivR1D1R1_Zero";
         case WORLANDINTG_R1_ZERO:
            return "R1_Zero";

         case WORLANDREDU:
            return "Reductor";
         case WORLANDREDU_ENERGY:
            return "Energy";
         case WORLANDREDU_ENERGYR2:
            return "Energy R2";
         case WORLANDREDU_ENERGYD1R1:
            return "Energy D1R1";
         case WORLANDREDU_POWER:
            return "Power";
         case WORLANDREDU_POWERR2:
            return "Power R2";
         case WORLANDREDU_POWERD1R1:
            return "Power D1R1";

         // Default output for unknown break point
         default:
            return "Unknown break point";
      }
   }

}
}
}
