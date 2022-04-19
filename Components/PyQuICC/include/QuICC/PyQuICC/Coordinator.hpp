/**
 * @file Coordinator.hpp
 * @brief Small coordinator for the Python interpreter to initialise and finalize with multiple wrappers
 */

#ifndef QUICC_PYQUICC_COORDINATOR_HPP
#define QUICC_PYQUICC_COORDINATOR_HPP

// First include
//
#include "QuICC/PyQuICC/SystemHeader.hpp"

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

namespace PyQuICC {

   /**
    * @brief Small coordinator for the Python interpreter to initialise and finalize with multiple wrappers
    */
   class Coordinator
   {
      public:
         /**
          * @brief Initialise the Python interpreter
          */
         static void init();

         /**
          * @brief Register Python wrapper
          */
         static void registerWrapper();

         /**
          * @brief Unregister Python wrapper
          */
         static void unregisterWrapper();

         /**
          * @brief Finalize the Python interpreter
          */
         static void finalize();
         
      protected:

      private:
         /**
          * @brief Counter for the number of active Python objects
          */
         static int sCounter; 

         /**
          * @brief Constructor
          */
         Coordinator();

         /**
          * @brief Destructor
          */
         ~Coordinator();
   };

}
}

#endif // QUICC_PYQUICC_COORDINATOR_HPP
