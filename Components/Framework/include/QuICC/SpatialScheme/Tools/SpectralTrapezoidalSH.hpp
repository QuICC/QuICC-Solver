/**
 * @file SpectralTrapezoidalSH.hpp
 * @brief Implementation of the tools for the trapezoidal + spherical harmonic spatial schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRAPEZOIDALSH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRAPEZOIDALSH_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBaseSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the trapezoidal + spherical harmonic spatial schemes
    */
   class SpectralTrapezoidalSH: public IBaseSH
   {
      public:
         /**
          * @brief Default ctor
          */
         SpectralTrapezoidalSH();

         /**
          * @brief ctor with explicit min truncation
          *
          * @param min  Minimal truncation
          */
         SpectralTrapezoidalSH(const int min);

         /**
          * @brief Default dtor
          */
         ~SpectralTrapezoidalSH() = default;

         /**
          * @brief Compute trapezoidal forward truncation
          */
         int truncationFwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute trapezoidal backward truncation
          */
         int truncationBwd(const int nN, const int j, const int k) final;

         /**
          * @brief Compute index
          */
         int index(const int nN, const int j, const int k) final;

         /**
          * @brief Check if chosen resolution is optimal
          */
         virtual bool isOptimal(const int nN, const int maxL);

      private:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_SPECTRALTRAPEZOIDALSH_HPP
