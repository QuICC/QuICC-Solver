/**
 * @file xLFlBuilder.hpp
 * @brief Implementation of the generic radial sphere + Spherical harmonics (Associated Legendre +  Fourier) scheme with spectral l ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_XLFLBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_XLFLBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHlBuilder.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the  generic radial + Spherical harmonics (Associated Legendre +  Fourier) scheme with spectral l ordering
    */
   class xLFlBuilder: public IRegularSHlBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit xLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~xLFlBuilder() = default;

         /**
          * @brief Add the transform setups to resolution
          */
         void addTransformSetups(SharedResolution spRes) const override;

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

#endif // QUICC_SPATIALSCHEME_3D_XLFLBUILDER_HPP
