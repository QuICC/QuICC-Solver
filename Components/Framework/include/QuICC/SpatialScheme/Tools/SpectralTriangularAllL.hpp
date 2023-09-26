/**
 * @file SpectralTriangularAllL.hpp
 * @brief Implementation of the tools for the triangular + spherical harmonics schemes with all harmonic degrees gathered
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARALLL_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARALLL_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/SpectralTrapezoidalAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes with all harmonic degrees gathered
    */
   class SpectralTriangularAllL: public SpectralTrapezoidalAllL
   {
      public:
         /**
          * @brief Default ctor
          */
         SpectralTriangularAllL() = default;

         /**
          * @brief ctor with explicit min truncation
          *
          * @param min  Minimal truncation
          */
         SpectralTriangularAllL(const int min);

         /**
          * @brief Default dtor
          */
         ~SpectralTriangularAllL() = default;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      private:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARALLL_HPP
