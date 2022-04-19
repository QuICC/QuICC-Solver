/**
 * @file StageTimer.hpp
 * @brief Implementation of a very simple stage timer 
 */

#ifndef QUICC_STAGETIMER_HPP
#define QUICC_STAGETIMER_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Timers/TimerMacro.h"

namespace QuICC {

   /**
    * @brief Implementation of a very simple stage timer
    */
   class StageTimer
   {
      public:
         /**
          * @brief Constructor
          */
         explicit StageTimer(const int digits = 1);

         /**
          * @brief Destructor
          */
         ~StageTimer();

         /**
          * @brief Set static IO flag
          */
         static void allowIo(const bool doesIo);

         /**
          * @brief New stage message
          *
          * @param msg Message to print
          */
         static void stage(const std::string& msg);

         /**
          * @brief Stage completed message
          *
          * @param msg Message to print
          */
         static void completed(const std::string& msg);

         /**
          * @brief Stage message
          *
          * @param msg Message to print
          */
         static void msg(const std::string& msg, const int space = 8);

         /**
          * @brief Start stage timing
          *
          * @param msg Message to print
          */
         void start(const std::string& msg, const int level = 0);

         /**
          * @brief End stage
          */
         void done();
         
      protected:

      private:
         /**
          * @brief Is allowed to do IO
          */
         static bool sDoesIo;

         /**
          * @brief 10^#digits to show
          */
         const double mcDigits;

         /**
          * @brief Level
          */
         int mLevel;

         /**
          * @brief The actual timer
          */
         TimerMacro mTimer;
   };

}

#endif // QUICC_STAGETIMER_HPP
