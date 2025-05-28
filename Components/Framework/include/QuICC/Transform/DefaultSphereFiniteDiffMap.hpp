/**
 * @file DefaultSphereFiniteDiffMap.hpp
 * @brief Default transform operator map in a sphere
 */

#ifndef QUICC_TRANSFORM_DEFAULTSPHEREFINITEDIFFMAP_HPP
#define QUICC_TRANSFORM_DEFAULTSPHEREFINITEDIFFMAP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/IOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Finite Differences transform in a sphere
    */
   class DefaultSphereFiniteDiffMap: public ITransformMap<FiniteDiff::Sphere::IOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultSphereFiniteDiffMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultSphereFiniteDiffMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTSPHEREFINITEDIFFMAP_HPP
