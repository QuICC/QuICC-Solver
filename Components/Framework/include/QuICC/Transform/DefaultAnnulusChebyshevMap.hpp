/**
 * @file DefaultAnnulusChebyshevMap.hpp
 * @brief Default transform operator map in a annulus
 */

#ifndef QUICC_TRANSFORM_DEFAULTANNULUSCHEBYSHEVMAP_HPP
#define QUICC_TRANSFORM_DEFAULTANNULUSCHEBYSHEVMAP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Chebyshev transform in a annulus
    */
   class DefaultAnnulusChebyshevMap: public ITransformMap<Fft::Chebyshev::IChebyshevOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultAnnulusChebyshevMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultAnnulusChebyshevMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTANNULUSCHEBYSHEVMAP_HPP
