/**
 * @file SphericalLorentzAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALLORENTZANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALLORENTZANELASTIC_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

#include "DenseSM/IGenericProfile.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalLorentzAnelastic
   {
      public:

         /**
          * @brief Set Lorentz term to S
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Lorentz term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Lorentz term to S
          */
         static void test(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Array& thGrid,    // Add theta grid
                         const Array& phGrid,    // Add phi grid 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalLorentzAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalLorentzAnelastic() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALLORENTZANELASTIC_HPP
