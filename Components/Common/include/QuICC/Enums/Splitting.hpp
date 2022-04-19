/**
 * @file Splitting.hpp
 * @brief Definition of some useful enums for splitting algorithms 
 */

#ifndef QUICC_SPLITTING_HPP
#define QUICC_SPLITTING_HPP

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
    * @brief Namespace holding the splitting algorithm related enums
    */
   namespace Splitting {

      /**
       * @brief Struct holding the possible splitting locations
       */
      struct Locations
      {
         /**
          * @name Enum for the splitting location IDs
          */
         enum Id {
            /// No splitting
            NONE,
            /// Splitting is done between first/second transform
            FIRST,
            /// Splitting is done between second/third transform
            SECOND,
            /// Splitting is done between first/second and second/third transform
            BOTH,
            /// Splitting is ONLY on slowest dimension for the first transform,
            COUPLED2D,
         };
      };

      /**
       * @brief Simple struct holding the parallelisation algorithms IDs
       */
      struct Algorithms {

         /**
          * @name Enum for algorithm IDs
          */
         enum Id {
            /// Serial code
            SERIAL,
            /// Single load splitting on first transform
            SINGLE1D,
            /// Single load splitting on second transform
            SINGLE2D,
            /// Load splitting on first and second transform (tubular, pencils)
            TUBULAR,
            /// Single load splitting on slowest direction on first transform
            COUPLED2D,
         };
      };

      /**
       * @brief Simple struct holding the transform grouper IDs
       */
      struct Groupers {

         /**
          * @name Enum for grouper IDs
          */
         enum Id {
            /// Group transforms by equation
            EQUATION,
            /// Group transforms at first transform stage
            SINGLE1D,
            /// Group transforms at second transform stage
            SINGLE2D,
            /// Group transforms at each transform stage
            TRANSFORM
         };
      };
   }
}

#endif // QUICC_SPLITTING_HPP
