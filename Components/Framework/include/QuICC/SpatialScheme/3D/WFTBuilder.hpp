/**
 * @file WFTBuilder.hpp
 * @brief Implementation of the cylindrical Worland(poly) + Fourier + Chebyshev(FFT) scheme
 */

#ifndef QUICC_SPATIALSCHEME_3D_WFTBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_WFTBUILDER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegular3DBuilder.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of cylindrical Worland(poly) + Fourier + Chebyshev(FFT) scheme
    */
   class WFTBuilder: public IRegular3DBuilder
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
          */
         explicit WFTBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~WFTBuilder();

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
         Transform::Poly::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup3D(SharedResolution spRes) const;
   };

}
}

#endif // QUICC_SPATIALSCHEME_3D_WFTBUILDER_HPP
