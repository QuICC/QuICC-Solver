/**
 * @file SphericalZComponent.hpp
 * @brief Implementation of the spherical Z component of a field
 */

#ifndef QUICC_PHYSICAL_SPHERICALZCOMPONENT_HPP
#define QUICC_PHYSICAL_SPHERICALZCOMPONENT_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical Z component of a field
    */
   class SphericalZComponent
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalZComponent() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalZComponent() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALZCOMPONENT_HPP
