/**
 * @file SphericalViscousDissipationAnelastic.hpp
 * @brief Implementation of the viscous dissipation term in the anelastic approximation
 * 
 * Q_\nu  
 */

#ifndef QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP

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
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of a generic primitive mass flux advection in 3D space
    */
   template <FieldComponents::Physical::Id TONE, 
             FieldComponents::Physical::Id TTWO, 
             FieldComponents::Physical::Id TTHREE> class SphericalViscousDissipationAnelastic
   {
      public:
         /**
          * @brief Set S to viscous dissipation
          */
          static void set(Framework::Selector::PhysicalScalarField &rS, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                          FieldComponents::Physical::Id> &u, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                          FieldComponents::Physical::Id> &gradQ, 
                          const MHDFloat c = 1.0);

         /**
          * @brief Add viscous dissipation to S
          */
          static void add(Framework::Selector::PhysicalScalarField &rS, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                          FieldComponents::Physical::Id> &u, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                          FieldComponents::Physical::Id> &gradQ, 
                          const MHDFloat c = 1.0);

         /**
          * @brief Substract viscous dissipation from S
          */
          static void sub(Framework::Selector::PhysicalScalarField &rS, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                          FieldComponents::Physical::Id> &u,
                           const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                           FieldComponents::Physical::Id> &gradQ, 
                           const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalViscousDissipationAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalViscousDissipationAnelastic() = default;
   };

   template <FieldComponents::Physical::Id TONE, 
             FieldComponents::Physical::Id TTWO, 
             FieldComponents::Physical::Id TTHREE> 
             void SphericalViscousDissipationAnelastic<TONE,TTWO,TTHREE>::set(Framework::Selector::PhysicalScalarField &rS, 
                                                                              const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, 
                                                                              const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, 
                                                                              const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.setData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, 
             FieldComponents::Physical::Id TTWO, 
             FieldComponents::Physical::Id TTHREE> 
             void SphericalViscousDissipationAnelastic<TONE,TTWO,TTHREE>::add(Framework::Selector::PhysicalScalarField &rS, 
                                                                              const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, 
                                                                              const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, 
                                                                              const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.addData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, 
             FieldComponents::Physical::Id TTWO, 
             FieldComponents::Physical::Id TTHREE> 
             void SphericalViscousDissipationAnelastic<TONE,TTWO,TTHREE>::sub(Framework::Selector::PhysicalScalarField &rS, 
                                                                              const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, 
                                                                              const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, 
                                                                              const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.subData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.subData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.subData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.subData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.subData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      }
   }

   }
}

#endif // QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP
