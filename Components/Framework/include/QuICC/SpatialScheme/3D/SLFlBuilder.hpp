/**
 * @file SLFlBuilder.hpp
 * @brief Implementation of the shell Chebyshev(FFT) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral l ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_SLFLBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_SLFLBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHlBuilder.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the shell Chebyshev(FFT) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral l ordering
    */
   class SLFlBuilder: public IRegularSHlBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          */
         explicit SLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~SLFlBuilder() = default;

         /**
          * @brief Add the transform setups to resolution
          */
         void addTransformSetups(SharedResolution spRes) const final;

      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         void setDimensions() final;

         /**
          * @brief Set transform costs
          */
         void setCosts() final;

         /**
          * @brief Set transform scalings
          */
         void setScalings() final;

         /**
          * @brief Set transform memory footprint
          */
         void setMemoryScore() final;

      private:
         /**
          * @brief Construct setup object for first transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Poly::ALegendre::SharedSetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup3D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_SLFLBUILDER_HPP
