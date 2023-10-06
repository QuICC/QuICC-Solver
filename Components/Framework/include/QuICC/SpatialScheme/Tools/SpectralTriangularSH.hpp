/**
 * @file SpectralTriangularSH.hpp
 * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARSH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARSH_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/SpectralTrapezoidalSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes
    */
   class SpectralTriangularSH: public SpectralTrapezoidalSH
   {
      public:
         /**
          * @brief Default ctor
          */
         SpectralTriangularSH() = default;

         /**
          * @brief ctor with explicit min truncation
          *
          * @param min Minimum truncation
          */
         SpectralTriangularSH(const int min);

         /**
          * @brief Default dtor
          */
         ~SpectralTriangularSH() = default;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      private:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARSH_HPP
