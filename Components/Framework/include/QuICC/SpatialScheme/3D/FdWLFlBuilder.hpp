/**
 * @file FdWLFlBuilder.hpp
 * @brief Implementation of the Finite Differences sphere + Spherical harmonics (Associated Legendre +  Fourier) scheme with spectral l ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_FDWLFLBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_FDWLFLBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/xLFlBuilder.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the Finite Difference sphere + Spherical harmonics (Associated Legendre +  Fourier) scheme with spectral l ordering
    */
   class FdWLFlBuilder: public xLFlBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          * @param options Scheme options
          */
         explicit FdWLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~FdWLFlBuilder() = default;

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

#endif // QUICC_SPATIALSCHEME_3D_FDWLFLBUILDER_HPP
