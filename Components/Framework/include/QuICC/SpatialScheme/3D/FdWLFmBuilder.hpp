/**
 * @file FdWLFmBuilder.hpp
 * @brief Implementation of the Finite Differences sphere + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_FDWLFMBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_FDWLFMBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/xLFmBuilder.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the Finite Differences sphere + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
    */
   class FdWLFmBuilder: public xLFmBuilder
   {
      public:
         /**
          * @brief Tune the shared resolution used by simulation
          */
         void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr) final;

         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit FdWLFmBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~FdWLFmBuilder() = default;

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
         Transform::SharedTransformSetup  spSetup1D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_FDWLFMBUILDER_HPP
