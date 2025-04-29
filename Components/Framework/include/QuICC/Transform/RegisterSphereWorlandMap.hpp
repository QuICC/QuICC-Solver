/**
 * @file RegisterSphereWorlandMap.hpp
 * @brief Register transform operators of the Worland transform in a sphere
 */

#ifndef QUICC_TRANSFORM_REGISTERSPHEREWORLANDMAP_HPP
#define QUICC_TRANSFORM_REGISTERSPHEREWORLANDMAP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Register transform operators of the Worland transform in a sphere
    */
   class RegisterSphereWorlandMap
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
         RegisterSphereWorlandMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterSphereWorlandMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERSPHEREWORLANDMAP_HPP
