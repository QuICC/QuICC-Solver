/**
 * @file TTTBuilder.hpp
 * @brief Implementation of the Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme
 */

#ifndef QUICC_SPATIALSCHEME_3D_TTTBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_TTTBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegular3DBuilder.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme
    */
   class TTTBuilder: public IRegular3DBuilder
   {
      public:
         /**
          * @brief Tune the shared resolution used by simulation
          */
         void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit TTTBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~TTTBuilder() = default;

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

         /**
          * @brief Construct setup object for third transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup3D(SharedResolution spRes) const;
   };

}
}

#endif // QUICC_SPATIALSCHEME_3D_TTTBUILDER_HPP
