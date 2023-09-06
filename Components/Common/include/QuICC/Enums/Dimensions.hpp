/**
 * @file Dimensions.hpp
 * @brief Definition of some useful enums of the dimensions of the model 
 */

#ifndef QUICC_DIMENSIONS_HPP
#define QUICC_DIMENSIONS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

namespace Dimensions {

   /**
    * @brief Simple struct holding the IDs of transform spaces
    */
   struct Transform {
      /**
       * @name Enums for the transform spaces
       */
      enum Id {
         /// First transform space
         TRA1D = 0,
         /// Second transform space
         TRA2D,
         /// Third transform space
         TRA3D,
         /// Spectral space
         SPECTRAL,
      };
   };

   /**
    * @brief Simple struct holding the IDs of simulation dimensions
    */
   struct Simulation {
      /**
       * @name Enums for the simulation dimensions
       */
      enum Id {
         /// First dimension data
         SIM1D = 0,
         /// Second dimension of data
         SIM2D,
         /// Third dimension data
         SIM3D
      };
   };

   /**
    * @brief Simple struct holding the IDs of data dimensions
    */
   struct Data {
      /**
       * @name Enums for the data dimensions
       */
      enum Id {
         /// First dimension of data for forward transform
         DATF1D = 0,
         /// First dimension of data for backward transform
         DATB1D,
         /// Second dimension of data
         DAT2D,
         /// Third dimension data
         DAT3D,
      };
   };

   /**
    * @brief Simple struct holding IDs for the different "spaces"
    */
   struct Space {
      /**
       * @name Enums for the dimension spaces
       */
      enum Id {
         /// Spectral space
         SPECTRAL = 0,
         /// Physical space
         PHYSICAL,
         /// Transform space (useful if spectral != transform)
         TRANSFORM
      };
   };
}
}

#endif // QUICC_DIMENSIONS_HPP
