/**
 * @file DefaultCartesianChebyshevMap.hpp
 * @brief Default transform operator map in a cartesian box
 */

#ifndef QUICC_TRANSFORM_DEFAULTCARTESIANCHEBYSHEVMAP_HPP
#define QUICC_TRANSFORM_DEFAULTCARTESIANCHEBYSHEVMAP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Chebyshev transform in a cartesian box
    */
   class DefaultCartesianChebyshevMap: public ITransformMap<Fft::Chebyshev::IChebyshevOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultCartesianChebyshevMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultCartesianChebyshevMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTCARTESIANCHEBYSHEVMAP_HPP
