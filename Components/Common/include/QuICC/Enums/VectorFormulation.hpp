/**
 * @file VectorFormulation.hpp
 * @brief Definition of some useful enums for the spatial formulation of vector fields 
 */

#ifndef QUICC_VECTORFORMULATION_HPP
#define QUICC_VECTORFORMULATION_HPP

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
    * @brief Formulation of spatial vector fields
    */
   struct VectorFormulation {

      /**
      * @name Enum for formulation of spatial vector fields
      */
      enum Id {
         /// Primitive formulation
         PRIMITIVE = 0,
         /// Toroidal/Poloidal formulation
         TORPOL,
         /// Spheroidal/Toroidal QST formulation
         QST,
      };
   };
}

#endif // QUICC_VECTORFORMULATION_HPP
