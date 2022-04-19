/**
 * @file ITimer.hpp
 * @brief Implementation of a timer interface 
 */

#ifndef QUICC_ITIMER_HPP
#define QUICC_ITIMER_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @brief Implementation of a timer interface
    */
   class ITimer
   {
      public:
         /**
          * @brief Constructor
          */
         ITimer() = default;

         /**
          * @brief Destructor
          */
         virtual ~ITimer() = default;

         /**
          * @brief Start clock
          */
         virtual void start() = 0;

         /**
          * @brief Stop clock
          */
         virtual void stop() = 0;

         /**
          * @brief Get elapsed time
          */
         virtual double time() const = 0;

         /**
          * @brief Get current time without stoping timer
          */
         virtual double queryTime() const = 0;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual double resetTimer() = 0;
         
      protected:

      private:
   };

}

#endif // QUICC_ITIMER_HPP
