/**
 * @file FFBuilder.hpp
 * @brief Implementation of the Fourier + Fourier scheme
 */

#ifndef QUICC_SPATIALSCHEME_FFBUILDER_HPP
#define QUICC_SPATIALSCHEME_FFBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/2D/IRegular2DBuilder.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Fourier + Fourier scheme
    */
   class FFBuilder: public IRegular2DBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit FFBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~FFBuilder() = default;

         /**
          * @brief Add the transform setups to resolution
          */
         void addTransformSetups(SharedResolution spRes) const final;

      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         void setDimensions() final;

      private:
         /**
          * @brief Construct setup object for first transform
          */
         Transform::Fft::Fourier::Complex::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup2D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_FFBUILDER_HPP
