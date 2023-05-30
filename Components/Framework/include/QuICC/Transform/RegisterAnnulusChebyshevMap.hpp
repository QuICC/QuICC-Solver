/**
 * @file RegisterAnnulusChebyshevMap.hpp
 * @brief Register transform operator map Chebyshev transform in a annulus
 */

#ifndef QUICC_TRANSFORM_REGISTERANNULUSCHEBYSHEVMAP_HPP
#define QUICC_TRANSFORM_REGISTERANNULUSCHEBYSHEVMAP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Register transform operator map Chebyshev transform in a annulus
    */
   class RegisterAnnulusChebyshevMap
   {
      public:
         typedef std::vector<std::shared_ptr<ITransformMap<Fft::Chebyshev::IChebyshevOperator> > > MapVector;

         /**
          * @brief Store transform operator to ID mapping
          */
         static MapVector& mapper();

      private:
         /**
          * @brief Constructor
          */
         RegisterAnnulusChebyshevMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterAnnulusChebyshevMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERANNULUSCHEBYSHEVMAP_HPP
