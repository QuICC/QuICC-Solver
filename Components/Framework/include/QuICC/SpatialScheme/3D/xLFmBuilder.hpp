/**
 * @file xLFmBuilder.hpp
 * @brief Implementation of the generic radial + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_XLFMBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_XLFMBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHmlBuilder.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the generic radial  + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
    */
   class xLFmBuilder: public IRegularSHmlBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit xLFmBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~xLFmBuilder() = default;

         /**
          * @brief Add the transform setups to resolution
          */
         void addTransformSetups(SharedResolution spRes) const override;

         /**
          * @brief Spectral ordering is different from transform
          */
         bool sameSpectralOrdering() const final;

      protected:

      private:
         /**
          * @brief Construct setup object for second transform
          */
         Transform::SharedTransformSetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::SharedTransformSetup  spSetup3D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_XLFMBUILDER_HPP
