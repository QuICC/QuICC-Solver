/**
 * @file ProfilerBase.hpp
 * @brief Profiling breakpoints
 */

#ifndef QUICC_DEBUG_PROFILER_BREAKPOINT_HPP
#define QUICC_DEBUG_PROFILER_BREAKPOINT_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Debug {

namespace Profiler {

   // Index for level 0 breakpoints
   static const unsigned int LVL0 = 1000000;
   // Index for level 1 breakpoints
   static const unsigned int LVL1 = LVL0/100;
   // Index for level 2 breakpoints
   static const unsigned int LVL2 = LVL1/100;
   // Index for level 3 breakpoints
   static const unsigned int LVL3 = LVL2/100;

   /**
    * @brief List of break points for the profiling
    */
   enum BreakPoint {
      // Black hole
      BLACKHOLE = 0,
      // Coarse profiling points (level 0)
      BWDTRANSFORM = 1*LVL0,
         BWDIN = BWDTRANSFORM + 1*LVL1,
         BWD1D = BWDTRANSFORM + 2*LVL1,
            BWD1DIN = BWD1D + 1*LVL2,
            BWD1DTRA = BWD1D + 2*LVL2,
               BWD1DTRA_A = BWD1DTRA + 1*LVL3,
               BWD1DTRA_B = BWD1DTRA + 2*LVL3,
               BWD1DTRA_C = BWD1DTRA + 3*LVL3,
               BWD1DTRA_D = BWD1DTRA + 4*LVL3,
               BWD1DTRA_E = BWD1DTRA + 5*LVL3,
               BWD1DTRA_F = BWD1DTRA + 6*LVL3,
               BWD1DTRA_G = BWD1DTRA + 7*LVL3,
               BWD1DTRA_H = BWD1DTRA + 8*LVL3,
               BWD1DTRA_I = BWD1DTRA + 9*LVL3,
            BWD1DOUT = BWD1D + 4*LVL2,
               BWD1DOUTWORK = BWD1DOUT + 1*LVL3,
               BWD1DOUTCOMM = BWD1DOUT + 2*LVL3,
         BWD2D = BWDTRANSFORM + 3*LVL1,
            BWD2DIN = BWD2D + 1*LVL2,
            BWD2DTRA = BWD2D + 2*LVL2,
            BWD2DOUT = BWD2D + 4*LVL2,
               BWD2DOUTWORK = BWD2DOUT + 1*LVL3,
               BWD2DOUTCOMM = BWD2DOUT + 2*LVL3,
         BWDND = BWDTRANSFORM + 4*LVL1,
            BWDNDIN = BWDND + 1*LVL2,
            BWDNDTRA = BWDND + 2*LVL2,
            BWDNDOUT = BWDND + 4*LVL2,
               BWDNDOUTWORK = BWDNDOUT + 1*LVL3,
               BWDNDOUTCOMM = BWDNDOUT + 2*LVL3,
         BWDOUT = BWDTRANSFORM + 5*LVL1,
      NONLINEAR = 2*LVL0,
      FWDTRANSFORM = 3*LVL0,
         FWDIN = FWDTRANSFORM + 1*LVL1,
         FWDND = FWDTRANSFORM + 2*LVL1,
            FWDNDIN = FWDND + 1*LVL2,
            FWDNDTRA = FWDND + 2*LVL2,
            FWDNDOUT = FWDND + 4*LVL2,
               FWDNDOUTWORK = FWDNDOUT + 1*LVL3,
               FWDNDOUTCOMM = FWDNDOUT + 2*LVL3,
         FWD2D = FWDTRANSFORM + 3*LVL1,
            FWD2DIN = FWD2D + 1*LVL2,
            FWD2DTRA = FWD2D + 2*LVL2,
            FWD2DOUT = FWD2D + 4*LVL2,
               FWD2DOUTWORK = FWD2DOUT + 1*LVL3,
               FWD2DOUTCOMM = FWD2DOUT + 2*LVL3,
         FWD1D = FWDTRANSFORM + 4*LVL1,
            FWD1DIN = FWD1D + 1*LVL2,
            FWD1DTRA = FWD1D + 2*LVL2,
            FWD1DOUT = FWD1D + 4*LVL2,
               FWD1DOUTWORK = FWD1DOUT + 1*LVL3,
               FWD1DOUTCOMM = FWD1DOUT + 2*LVL3,
         FWDOUT = FWDTRANSFORM + 5*LVL1,
      PROGNOSTICEQUATION = 4*LVL0,
         TSTEPEXIN = PROGNOSTICEQUATION + 1*LVL1,
         TSTEPIN = PROGNOSTICEQUATION + 2*LVL1,
         TSTEPSOLVE = PROGNOSTICEQUATION + 3*LVL1,
         TSTEPOUT = PROGNOSTICEQUATION + 4*LVL1,
      DIAGNOSTICEQUATION = 5*LVL0,
         DIAGNOSTICEXIN = DIAGNOSTICEQUATION + 1*LVL1,
         DIAGNOSTICSOLVE = DIAGNOSTICEQUATION + 2*LVL1,
      TRIVIALEQUATION = 6*LVL0,
         TRIVIALEXIN = TRIVIALEQUATION + 1*LVL1,
         TRIVIALSOLVE = TRIVIALEQUATION + 2*LVL1,
      CONTROL = 7*LVL0,
      IO = 8*LVL0,
      CONSTRAIN = 9*LVL0,
      // Details for transforms
      WORLANDTRA = 50*LVL0,
         WORLANDPROJ = WORLANDTRA + 1*LVL1,
            WORLANDPROJ_P = WORLANDPROJ + 1*LVL2,
            WORLANDPROJ_D1 = WORLANDPROJ + 2*LVL2,
            WORLANDPROJ_D1R1 = WORLANDPROJ + 3*LVL2,
            WORLANDPROJ_DIVR1D1R1 = WORLANDPROJ + 4*LVL2,
            WORLANDPROJ_DIVR1D1R1_ZERO = WORLANDPROJ + 5*LVL2,
            WORLANDPROJ_DIVR1 = WORLANDPROJ + 6*LVL2,
            WORLANDPROJ_DIVR1_ZERO = WORLANDPROJ + 7*LVL2,
            WORLANDPROJ_SPHLAPL = WORLANDPROJ + 8*LVL2,
         WORLANDINTG = WORLANDTRA + 2*LVL1,
            WORLANDINTG_P = WORLANDINTG + 1*LVL2,
            WORLANDINTG_R1 = WORLANDINTG + 2*LVL2,
            WORLANDINTG_DIVR1 = WORLANDINTG + 3*LVL2,
            WORLANDINTG_DIVR1D1R1 = WORLANDINTG + 4*LVL2,
            WORLANDINTG_I2 = WORLANDINTG + 5*LVL2,
            WORLANDINTG_I2DIVR1D1R1_ZERO = WORLANDINTG + 6*LVL2,
            WORLANDINTG_I2DIVR1_ZERO = WORLANDINTG + 7*LVL2,
            WORLANDINTG_I4DIVR1D1R1_ZERO = WORLANDINTG + 8*LVL2,
            WORLANDINTG_I4DIVR1_ZERO = WORLANDINTG + 9*LVL2,
            WORLANDINTG_P_ZERO = WORLANDINTG + 10*LVL2,
            WORLANDINTG_I2_ZERO = WORLANDINTG + 11*LVL2,
            WORLANDINTG_DIVR1_ZERO = WORLANDINTG + 12*LVL2,
            WORLANDINTG_DIVR1D1R1_ZERO = WORLANDINTG + 13*LVL2,
            WORLANDINTG_R1_ZERO = WORLANDINTG + 14*LVL2,
         WORLANDREDU = WORLANDTRA + 3*LVL1,
            WORLANDREDU_ENERGY = WORLANDREDU + 1*LVL2,
            WORLANDREDU_ENERGYR2 = WORLANDREDU + 2*LVL2,
            WORLANDREDU_ENERGYD1R1 = WORLANDREDU + 3*LVL2,
            WORLANDREDU_POWER = WORLANDREDU + 4*LVL2,
            WORLANDREDU_POWERR2 = WORLANDREDU + 5*LVL2,
            WORLANDREDU_POWERD1R1 = WORLANDREDU + 6*LVL2,
   };

   /**
    * @brief Shift breakpoint by integer
    */
   BreakPoint shiftPoint(BreakPoint point, const int shift);

   /**
    * @brief Get a human readable name for break point
    *
    * @param point Break point
    */
   std::string pointName(BreakPoint point);

}
}
}

#endif // QUICC_DEBUG_PROFILER_BREAKPOINT_HPP
