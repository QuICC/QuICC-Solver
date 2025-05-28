/**
 * @file RegisterSphereFiniteDiffMap.hpp
 * @brief Register transform operators of the Finite Differnces transform in a sphere
 */

#ifndef QUICC_TRANSFORM_REGISTERSPHEREFINITEDIFFMAP_HPP
#define QUICC_TRANSFORM_REGISTERSPHEREFINITEDIFFMAP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/IOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Register transform operators of the Finite Differences transform in a sphere
    */
   class RegisterSphereFiniteDiffMap
   {
      public:
         typedef std::vector<std::shared_ptr<ITransformMap<FiniteDiff::Sphere::IOperator> > > MapVector;

         /**
          * @brief Store transform operator to ID mapping
          */
         static MapVector& mapper();

      private:
         /**
          * @brief Constructor
          */
         RegisterSphereFiniteDiffMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterSphereFiniteDiffMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERSPHEREFINITEDIFFMAP_HPP
