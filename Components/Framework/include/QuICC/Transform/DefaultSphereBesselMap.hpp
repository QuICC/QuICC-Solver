/**
 * @file DefaultSphereBesselMap.hpp
 * @brief Default transform operator map for spherical Bessel expansion in a sphere
 */

#ifndef QUICC_TRANSFORM_DEFAULTSPHEREBESSELMAP_HPP
#define QUICC_TRANSFORM_DEFAULTSPHEREBESSELMAP_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/ITransformMap.hpp"
#include "QuICC/Transform/Poly/Bessel/IBesselOperator.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Default transform operator map for spherical Bessel expansion in a sphere
    */
   class DefaultSphereBesselMap: public ITransformMap<Poly::Bessel::IBesselOperator>
   {
      public:
         /**
          * @brief Constructor
          */
         DefaultSphereBesselMap() = default;

         /**
          * @brief Destructor
          */
         virtual ~DefaultSphereBesselMap() = default;

         /**
          * @brief Store transform operator to ID mapping
          */
         void operator()(MapType& m) const override;

      private:
   };

} // Transform
} // QuICC

#endif // QUICC_TRANSFORM_DEFAULTSPHEREBESSELMAP_HPP
