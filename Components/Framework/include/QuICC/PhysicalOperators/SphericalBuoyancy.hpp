/**
 * @file SphericalBuoyancy.hpp
 * @brief Implementation of the spherical buoyancy term
 */

#ifndef QUICC_PHYSICAL_SPHERICALBUOYANCY_HPP
#define QUICC_PHYSICAL_SPHERICALBUOYANCY_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical buoyancy term
    */
   class SphericalBuoyancy
   {
      public:
         /**
          * @brief Set S 
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Framework::Selector::PhysicalScalarField &q, const MHDFloat c = 1.0);

         /**
          * @brief Add to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Framework::Selector::PhysicalScalarField &q, const MHDFloat c = 1.0);

         /**
          * @brief Substract from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Framework::Selector::PhysicalScalarField &v, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalBuoyancy() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalBuoyancy() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALBUOYANCY_HPP
