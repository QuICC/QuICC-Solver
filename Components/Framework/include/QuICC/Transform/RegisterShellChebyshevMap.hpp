/**
 * @file RegisterShellChebyshevMap.hpp
 * @brief Register transform operators of the Chebyshev transform in a spherical shell
 */

#ifndef QUICC_TRANSFORM_REGISTERSHELLCHEBYSHEVMAP_HPP
#define QUICC_TRANSFORM_REGISTERSHELLCHEBYSHEVMAP_HPP

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
    * @brief Register transform operators of the Chebyshev transform in a spherical shell
    */
   class RegisterShellChebyshevMap
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
         RegisterShellChebyshevMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterShellChebyshevMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERSHELLCHEBYSHEVMAP_HPP
