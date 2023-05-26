/**
 * @file WLFlBuilder.hpp
 * @brief Implementation of the sphere Worland + Spherical harmonics (Associated Legendre +  Fourier) scheme with spectral l ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_WLFLBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_WLFLBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHlBuilder.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the sphere Worland + Spherical harmonics (Associated Legendre +  Fourier) scheme with spectral l ordering
    */
   class WLFlBuilder: public IRegularSHlBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          */
         explicit WLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~WLFlBuilder() = default;

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
         Transform::SharedTransformSetup  spSetup1D(SharedResolution spRes) const;

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

#endif // QUICC_SPATIALSCHEME_3D_WLFLBUILDER_HPP
