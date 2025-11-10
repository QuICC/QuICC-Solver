/**
 * @file SphericalCoriolisAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

#include "DenseSM/IGenericProfile.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalCoriolisAnelastic
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Array& cosTheta, 
                         const Array& sinTheta, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, 
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Array& cosTheta, 
                         const Array& sinTheta, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, 
                         const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Array& cosTheta, 
                         const Array& sinTheta, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, 
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalCoriolisAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalCoriolisAnelastic() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP
