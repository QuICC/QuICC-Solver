/**
 * @file TTBuilder.hpp
 * @brief Implementation of the Chebyshev(FFT) + Chebyshev(FFT) scheme
 */

#ifndef QUICC_SPATIALSCHEME_TTBUILDER_HPP
#define QUICC_SPATIALSCHEME_TTBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/2D/IRegular2DBuilder.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Chebyshev(FFT) + Chebyshev(FFT) scheme
    */
   class TTBuilder: public IRegular2DBuilder
   {
      public:
         /**
          * @brief Tune the shared resolution used by simulation
          */
         void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr) override;

         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit TTBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~TTBuilder() = default;

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
         Transform::Fft::Chebyshev::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup2D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TTBUILDER_HPP
