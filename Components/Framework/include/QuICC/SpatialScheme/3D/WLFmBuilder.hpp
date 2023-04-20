/**
 * @file WLFmBuilder.hpp
 * @brief Implementation of the sphere Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_WLFMBUILDER_HPP
#define QUICC_SPATIALSCHEME_3D_WLFMBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHmlBuilder.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the sphere Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
    */
   class WLFmBuilder: public IRegularSHmlBuilder
   {
      public:
         /**
          * @brief Tune the shared resolution used by simulation
          */
         void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr) override;

         /**
          * @brief Constructor
          *
          * @param dim     Spectral dimensions
          * @param purpose Grid purpose
          */
         explicit WLFmBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~WLFmBuilder() = default;

         /**
          * @brief Add the transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const override;

         /**
          * @brief Spectral ordering is different from transform
          */
         bool sameSpectralOrdering() const override;

      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions() override;

         /**
          * @brief Set transform costs
          */
         virtual void setCosts() override;

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings() override;

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemoryScore() override;

      private:
         /**
          * @brief Construct setup object for first transform
          */
         Transform::Poly::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Poly::ALegendre::SharedSetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup3D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_WLFMBUILDER_HPP
