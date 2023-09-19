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
#include "QuICC/SpatialScheme/Tools/IBaseAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the triangular + spherical harmonic spatial schemes with all harmonic degrees gathered
    */
   class SpectralTriangularAllL: public IBaseAllL
   {
      public:
         /**
          * @brief Default ctor
          */
         SpectralTriangularAllL() = default;

         /**
          * @brief Default dtor
          */
         ~SpectralTriangularAllL() = default;

         /**
          * @brief Compute backward truncation
          */
         int truncationBwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute index
          */
         int index(const int nN, const int j, const int k) final;

         /**
          * @brief Check if chosen resolution is optimal
          */
         bool isOptimal(const int nN, const int maxL) final;

      private:
         /**
          * @brief Minimal truncation for highest modes
          */
         static constexpr const int MIN_TRUNCATION = 3;
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRIANGULARALLL_HPP
