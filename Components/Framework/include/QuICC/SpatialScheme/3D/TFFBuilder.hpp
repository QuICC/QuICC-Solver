/**
 * @file TFFBuilder.hpp
 * @brief Implementation of the Chebyshev(FFT) + Fourier + Fourier scheme
 */

#ifndef QUICC_SPATIALSCHEME_3D_TFFBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_TFFBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegular3DBuilder.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Chebyshev(FFT) + Fourier + Fourier scheme
    */
   class TFFBuilder: public IRegular3DBuilder
   {
      public:
         /**
          * @brief Interpret the configuration dimensions
          */
         static void interpretConfigDimensions(ArrayI& rDim);

         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit TFFBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~TFFBuilder() = default;

         /**
          * @brief Add the transform setups to resolution
          */
         void addTransformSetups(SharedResolution spRes) const final;

         /**
          * @brief Add index counter to shared resolution
          */
         void addIndexCounter(SharedResolution spRes) final;

      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         void setDimensions() final;

      private:
         /**
          * @brief Construct setup object for first transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Fft::Fourier::Complex::SharedSetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup3D(SharedResolution spRes) const;

   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_TFFBUILDER_HPP
