/**
 * @file RegisterCylinderWorlandMap.hpp
 * @brief Register transform operators of the Worland transform in a cylinder
 */

#ifndef QUICC_TRANSFORM_REGISTERCYLINDERWORLANDMAP_HPP
#define QUICC_TRANSFORM_REGISTERCYLINDERWORLANDMAP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Register transform operators of the Worland transform in a cylinder
    */
   class RegisterCylinderWorlandMap
   {
      public:
         typedef std::vector<std::shared_ptr<ITransformMap<Poly::Worland::IWorlandOperator> > > MapVector;

         /**
          * @brief Store transform operator to ID mapping
          */
         static MapVector& mapper();

      private:
         /**
          * @brief Constructor
          */
         RegisterCylinderWorlandMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterCylinderWorlandMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERCYLINDERWORLANDMAP_HPP
