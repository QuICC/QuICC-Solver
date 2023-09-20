/**
 * @file DefaultCylinderWorlandMap.hpp
 * @brief Default transform operator map in a cylinder
 */

#ifndef QUICC_TRANSFORM_DEFAULTCYLINDERWORLANDMAP_HPP
#define QUICC_TRANSFORM_DEFAULTCYLINDERWORLANDMAP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Worland transform in a cylinder
    */
   class DefaultCylinderWorlandMap: public ITransformMap<Poly::Worland::IWorlandOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultCylinderWorlandMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultCylinderWorlandMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTCYLINDERWORLANDMAP_HPP
