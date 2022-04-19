/**
 * @file FieldIds.hpp
 * @brief Definition of some useful enums used to access fields by ID 
 */

#ifndef QUICC_FIELDIDS_HPP
#define QUICC_FIELDIDS_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @brief Simple struct to hold the field types
    */
   struct FieldType
   {
      /**
       * @brief Enum for the field types
       */
      enum Id {
         /// Scalar field
         SCALAR,
         /// Vector field
         VECTOR,
         /// Gradient field
         GRADIENT,
         /// Curl field
         CURL,
         /// 2nd order gradient field
         GRADIENT2,
         /// divergence field
         DIVERGENCE,
      };
   };

   /**
    * @brief Simple struct to hold the field components
    */
   struct FieldComponents
   {
      /**
       * @brief Struct for the physical field components
       */
      struct Physical
      {
         /**
          * @brief Enum for physical field vector components
          */
         enum Id {
            /// X component of cartesian field
            X,
            /// Radial component of cylindrical or spherical field
            R,
            /// Y component of cartesian field
            Y,
            /// Theta component of cylindrical or spherical field
            THETA,
            /// Z component of cartesian field
            Z,
            /// Phi component of spherical field
            PHI,

            /// Is a scalar
            SCALAR,

            /// Is not used
            NOTUSED,
         };
      };

      /**
       * @brief Struct for the spectral field components
       */
      struct Spectral
      {
         /**
          * @brief Enum for Spectral field vector components
          */
         enum Id {
            /// X component of cartesian field
            X,
            /// Radial component of cylindrical or spherical field
            R,
            /// Toroidal component
            TOR,
            /// Q component of Spheroidal/Toroidal
            Q,
            /// Y component of cartesian field
            Y,
            /// Theta component of cylindrical or spherical field
            THETA,
            /// Poloidal component of Toroidal/Poloidal
            POL,
            /// S component of Spheroidal/Toroidal
            S,
            /// Z component of cartesian field
            Z,
            /// Phi component of spherical field
            PHI,
            /// T component of spheroidal/toroidal field
            T,

            /// Is spectral scalar
            SCALAR,

            /// Is not used
            NOTUSED,
         };
      };
   };

   /// Typedef for a full ID for a spectral field component
   typedef std::pair<std::size_t,FieldComponents::Spectral::Id>   SpectralFieldId;

   /// Typedef for a full ID for a physical field component
   typedef std::pair<std::size_t,FieldComponents::Physical::Id>   PhysicalFieldId;
}

#endif // QUICC_FIELDIDS_HPP
