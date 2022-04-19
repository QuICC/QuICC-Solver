/**
 * @file PerformanceTimer.hpp
 * @brief Implementation of a very simple stage timer with focus on performance
 */

#ifndef QUICC_PERFORMANCETIMER_HPP
#define QUICC_PERFORMANCETIMER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "Timers/TimerMacro.h"

namespace QuICC {

   /**
    * @brief Implementation of a very simple stage timer with focus on performance
    */
   class PerformanceTimer
   {
      public:
         /**
          * @brief Constructor
          */
         explicit PerformanceTimer(const int nTimer, const std::string& name = "", const int digits = 6);

         /**
          * @brief Constructor
          */
         explicit PerformanceTimer(const std::map<int,std::string>& tag, const std::string& name = "", const int digits = 6);

         /**
          * @brief Destructor
          */
         ~PerformanceTimer();

         /**
          * @brief Set static IO flag
          */
         static void allowIo(const bool doesIo);

         /**
          * @brief Initialize timers with ID 0..nTimer-1
          */
         void init(const int nTimer);

         /**
          * @brief Set name of timer
          */
         void setName(const std::string& name);

         /**
          * @brief Add timer
          */
         void add(const int id, const std::string& tag = "");

         /**
          * @brief Set tag of timer
          */
         void setTag(const int id, const std::string& tag);

         /**
          * @brief Start stage timing
          *
          * @param id ID of stage
          */
         void start(const int id);

         /**
          * @brief End stage
          *
          * @param id ID of stage
          */
         void stop(const int id);

         /**
          * @brief Get stage time
          */
         double time(const int id);

         /**
          * @brief Show all timings
          */
         void show();
         
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
          * @brief Show values upon destruction
          */
         bool mShowOnDelete;

         /**
          * @brief Name of timer
          */
         std::string mName;

         /**
          * @brief The actual timers
          */
         std::map<int,TimerMacro> mTimer;

         /**
          * @brief Tags for timer
          */
         std::map<int,std::string> mTag;
   };

}

#endif // QUICC_PERFORMANCETIMER_HPP
