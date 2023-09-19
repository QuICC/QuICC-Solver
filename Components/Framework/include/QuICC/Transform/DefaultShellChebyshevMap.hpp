/**
 * @file DefaultShellChebyshevMap.hpp
 * @brief Default transform operator map in a spherical shell
 */

#ifndef QUICC_TRANSFORM_DEFAULTSHELLCHEBYSHEVMAP_HPP
#define QUICC_TRANSFORM_DEFAULTSHELLCHEBYSHEVMAP_HPP

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
    * @brief Implementation of the Chebyshev transform in a spherical shell
    */
   class DefaultShellChebyshevMap: public ITransformMap<Fft::Chebyshev::IChebyshevOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultShellChebyshevMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultShellChebyshevMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTSHELLCHEBYSHEVMAP_HPP
