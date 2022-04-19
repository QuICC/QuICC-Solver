/**
 * @file GridPurpose.hpp
 * @brief Definition of some useful enums for the purpose of the grid values 
 */

#ifndef QUICC_GRIDPURPOSE_HPP
#define QUICC_GRIDPURPOSE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @brief Purpose of the grid values
    */
   struct GridPurpose {

      /**
      * @name Enum for purpose of the grid values
      */
      enum Id {
         /// Used (mainly) for nonlinear simulation
         SIMULATION = 0,
         /// Used (mainly) for visualization
         VISUALIZATION,
      };
   };
}

#endif // QUICC_GRIDPURPOSE_HPP
