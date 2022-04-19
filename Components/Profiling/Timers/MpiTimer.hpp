/**
 * @file MpiTimer.hpp
 * @brief Implementation of a MPI timer 
 */

#ifndef QUICC_MPITIMER_HPP
#define QUICC_MPITIMER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "ITimer.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a MPI timer
    */
   class MpiTimer: public ITimer
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         explicit MpiTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~MpiTimer();

         /**
          * @brief Start clock
          */
         virtual void start();

         /**
          * @brief Stop clock
          */
         virtual void stop();

         /**
          * @brief Get elapsed time
          */
         virtual double time() const;

         /**
          * @brief Query current time without changes to state
          */
         virtual double queryTime() const;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual double resetTimer();
         
      protected:

      private:
         /**
          * @brief Start
          */
         double mStart;

         /**
          * @brief Stop
          */
         double mStop;
   };

}

#endif // QUICC_MPITIMER_HPP
