/**
 * @file TimerMacro.h
 * @brief Preprocessor macros to setup timers depending on parallelisation options 
 */

#ifndef QUICC_TIMERMACRO_H
#define QUICC_TIMERMACRO_H

#ifdef QUICC_MPI
   // include MPI timer
   #include "MpiTimer.hpp"

   namespace QuICC {
      /// Typedef for a generic Timer based on the MpiTimer
      typedef MpiTimer  TimerMacro;
   }
#else
   // include serial timer
   #include "SerialTimer.hpp"

   namespace QuICC {
      /// Typedef for a generic Timer based on the SerialTimer
      typedef SerialTimer  TimerMacro;
   }
#endif // QUICC_MPI

#endif // QUICC_TIMERMACRO_H
