/**
 * @file StatisticCoordinator.hpp
 * @brief Coordinator for the statistics computations 
 */

#ifndef QUICC_STATISTICS_STATISTICCOORDINATOR_HPP
#define QUICC_STATISTICS_STATISTICCOORDINATOR_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Statistics {

   /**
    * @brief Coordinator for the statistics computations
    */
   class StatisticCoordinator
   {
      public:
         /**
          * @brief Constructor
          */
         StatisticCoordinator();

         /**
          * @brief Constructor
          */
         ~StatisticCoordinator();

         /**
          * @brief Initialise the coordinator
          */
         void init(); 

      protected:

      private:
   };
}
}

#endif // QUICC_STATISTICS_STATISTICCOORDINATOR_HPP
