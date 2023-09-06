/** 
 * @file FBuilder.hpp
 * @brief Implementation of the Fourier scheme
 */

#ifndef QUICC_SPATIALSCHEME_FBUILDER_HPP
#define QUICC_SPATIALSCHEME_FBUILDER_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/1D/IRegular1DBuilder.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Fourier scheme
    */
   class FBuilder: public IRegular1DBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations 
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit FBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         ~FBuilder() = default;

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
          * @brief Construct setup object for X transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup1D(SharedResolution spRes) const;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_FBUILDER_HPP
