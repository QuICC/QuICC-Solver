/** 
 * @file TTBuilder.hpp
 * @brief Implementation of the Chebyshev(FFT) + Chebyshev(FFT) scheme
 */

#ifndef QUICC_SPATIALSCHEME_TTBUILDER_HPP
#define QUICC_SPATIALSCHEME_TTBUILDER_HPP

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
#include "QuICC/SpatialScheme/2D/IRegular2DBuilder.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of Chebyshev(FFT) + Chebyshev(FFT) scheme
    */
   class TTBuilder: public IRegular2DBuilder
   {
      public:
         /**
          * @brief Tune the shared resolution used by simulation
          */
         void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Constructor
          *
          * @param dim Chebyshev truncations
          */
         explicit TTBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~TTBuilder();

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
          */
         virtual void setCosts();

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings();

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemoryScore();

      private:
         /**
          * @brief Construct setup object for first transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::Fft::Chebyshev::SharedSetup  spSetup2D(SharedResolution spRes) const;
   };

}
}

#endif // QUICC_SPATIALSCHEME_TTBUILDER_HPP
