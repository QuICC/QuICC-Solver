/**
 * @file DefaultSphereWorlandMap.hpp
 * @brief Default transform operator map in a sphere
 */

#ifndef QUICC_TRANSFORM_DEFAULTSPHEREWORLANDMAP_HPP
#define QUICC_TRANSFORM_DEFAULTSPHEREWORLANDMAP_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Worland transform in a sphere
    */
   class DefaultSphereWorlandMap: public ITransformMap<Poly::Worland::IWorlandOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultSphereWorlandMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultSphereWorlandMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTSPHEREWORLANDMAP_HPP
