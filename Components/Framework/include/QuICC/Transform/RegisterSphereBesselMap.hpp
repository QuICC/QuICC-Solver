/**
 * @file RegisterSphereBesselMap.hpp
 * @brief Register transform operators of the spherical Bessel transform in a sphere
 */

#ifndef QUICC_TRANSFORM_REGISTERSPHEREBESSELMAP_HPP
#define QUICC_TRANSFORM_REGISTERSPHEREBESSELMAP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Poly/Bessel/IBesselOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Register transform operators of the spherical Bessel transform in a sphere
    */
   class RegisterSphereBesselMap
   {
      public:
         typedef std::vector<std::shared_ptr<ITransformMap<Poly::Bessel::IBesselOperator> > > MapVector;

         /**
          * @brief Store transform operator to ID mapping
          */
         static MapVector& mapper();

      private:
         /**
          * @brief Constructor
          */
         RegisterSphereBesselMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~RegisterSphereBesselMap() = default;

   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_REGISTERSPHEREBESSELMAP_HPP
