/** 
 * @file SpectralUniformAllL.hpp
 * @brief Implementation of the tools for the uniform + spherical harmonics schemes with all harmonic degrees gathered
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SPECTRALUNIFORMALLL_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SPECTRALUNIFORMALLL_HPP

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
    * @brief Implementation of the tools for the uniform + spherical harmonic spatial schemes with all harmonic degrees gathered
    */
   class SpectralUniformAllL: public IBaseAllL
   {
      public:
         /**
          * @brief Default ctor
          */
         SpectralUniformAllL() = default;

         /**
          * @brief Default dtor
          */
         ~SpectralUniformAllL() = default;

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

      protected:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SPECTRALUNIFORMALLL_HPP
