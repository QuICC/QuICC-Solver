/**
 * @file TriangularSH.hpp
 * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_TRIANGULARSH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_TRIANGULARSH_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/TrapezoidalSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes
    */
   class TriangularSH: public TrapezoidalSH
   {
      public:
         /**
          * @brief Default ctor
          */
         TriangularSH() = default;

         /**
          * @brief ctor with explicit min truncation
          */
         TriangularSH(const int min);

         /**
          * @brief Default dtor
          */
         ~TriangularSH() = default;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      private:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_TRIANGULARSH_HPP
