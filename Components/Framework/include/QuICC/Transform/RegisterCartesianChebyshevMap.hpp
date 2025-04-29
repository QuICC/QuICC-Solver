/**
 * @file RegisterCartesianChebyshevMap.hpp
 * @brief Register transform operators of the Chebyshev transform in a cartesian box
 */

#ifndef QUICC_TRANSFORM_REGISTERCARTESIANCHEBYSHEVMAP_HPP
#define QUICC_TRANSFORM_REGISTERCARTESIANCHEBYSHEVMAP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Register transform operators of the Chebyshev transform in a cartesian box
    */
   class RegisterCartesianChebyshevMap
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
         RegisterCartesianChebyshevMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterCartesianChebyshevMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERCARTESIANCHEBYSHEVMAP_HPP
