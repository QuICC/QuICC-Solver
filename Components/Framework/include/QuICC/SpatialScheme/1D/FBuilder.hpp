/** 
 * @file FBuilder.hpp
 * @brief Implementation of the Fourier scheme
 */

#ifndef QUICC_SPATIALSCHEME_FBUILDER_HPP
#define QUICC_SPATIALSCHEME_FBUILDER_HPP

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
          * @param dim  Chebyshev truncation
          */
         explicit FBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~FBuilder();

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const;

         /**
          * @brief Add the transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const;
         
      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions();

         /**
          * @brief Set transform costs
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setCosts();

         /**
          * @brief Set transform scalings
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setScalings();

         /**
          * @brief Set transform memory footprint
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setMemoryScore();

      private:
         /**
          * @brief Construct setup object for X transform
          */
         Transform::Fft::Fourier::Mixed::SharedSetup  spSetup1D(SharedResolution spRes) const;
   };

}
}

#endif // QUICC_SPATIALSCHEME_FBUILDER_HPP
