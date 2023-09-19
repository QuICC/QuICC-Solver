/**
 * @file IRegularSHlBuilder.hpp
 * @brief Implementation of the Regular basis + Spherical Harmonics scheme with l spectral ordering
 */

#ifndef QUICC_SPATIALSCHEME_IREGULARSHLBUILDER_HPP
#define QUICC_SPATIALSCHEME_IREGULARSHLBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/I3DBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Regular basis + Spherical Harmonics scheme with l spectral ordering
    */
   class IRegularSHlBuilder: public I3DBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim     0: radial, 1: latitudinal, 2: longitudinal
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit IRegularSHlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         virtual ~IRegularSHlBuilder() = default;

         /**
          * @brief Get spatial scheme dimensions
          */
         virtual ArrayI resolution() const override;

         /**
          * @brief Add index counter to shared resolution
          */
         virtual void addIndexCounter(SharedResolution spRes) override;

      protected:
         /**
          * @brief Get truncation tool
          */
         virtual std::shared_ptr<Tools::IBase> truncationTools(const Dimensions::Transform::Id transId) const override;

         /**
          * @brief Regular truncation
          */
         int   mI;

         /**
          * @brief Spherical harmonic degree
          */
         int   mL;

         /**
          * @brief Spherical harmonic order
          */
         int   mM;

      private:
   };
}
}

#endif // QUICC_SPATIALSCHEME_IREGULARSHLBUILDER_HPP
