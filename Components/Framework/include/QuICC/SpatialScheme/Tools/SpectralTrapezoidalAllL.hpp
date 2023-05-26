/** 
 * @file SpectralTrapezoidalAllL.hpp
 * @brief Implementation of the tools for the trapezoidal + spherical harmonics schemes with all harmonic degrees gathered
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRAPEZOIDALALLL_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRAPEZOIDALALLL_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBaseAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the trapezoidal + spherical harmonic spatial schemes with all harmonic degrees gathered
    */
   class SpectralTrapezoidalAllL: public IBaseAllL
   {
      public:
         /**
          * @brief Default ctor
          */
         SpectralTrapezoidalAllL() = default;

         /**
          * @brief Default dtor
          */
         ~SpectralTrapezoidalAllL() = default;

         /**
          * @brief Compute backward truncation
          */
         int truncationBwd(const int nN, const int l) final;

         /**
          * @brief Compute index
          */
         int index(const int nN, const int k) final;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      private:
         /**
          * @brief Minimal truncation for highest modes
          */
         static const int MIN_TRUNCATION;
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRAPEZOIDALALLL_HPP
