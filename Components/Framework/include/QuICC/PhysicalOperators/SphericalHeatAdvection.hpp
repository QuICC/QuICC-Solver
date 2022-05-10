/**
 * @file SphericalHeatAdvection.hpp
 * @brief Implementation of a spherical internal heat advection
 */

#ifndef QUICC_PHYSICAL_SPHERICALHEATADVECTION_HPP
#define QUICC_PHYSICAL_SPHERICALHEATADVECTION_HPP

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
#include "QuICC/Framework/Selector/ScalarField.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of a spherical internal heating advection
    */
   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> class SphericalHeatAdvection
   {
      public:
         /**
          * @brief Set S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q+q_b)\f$
          */
          static void set(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q + q_b)\f$
          */
          static void add(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q + q_b)\f$
          */
          static void sub(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalHeatAdvection() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalHeatAdvection() = default;
   };

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void SphericalHeatAdvection<TONE,TTWO,TTHREE>::set(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      if(c != 1.0)
      {
         rS.setData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            rS.subSlice((c*r(iR_))*u.comp(TONE).slice(iR), iR);
         }

         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.setData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            rS.subSlice(r(iR_)*u.comp(TONE).slice(iR), iR);
         }

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void SphericalHeatAdvection<TONE,TTWO,TTHREE>::add(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      if(c != 1.0)
      {
         rS.addData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            rS.subSlice((c*r(iR_))*u.comp(TONE).slice(iR), iR);
         }

         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.addData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            rS.subSlice(r(iR_)*u.comp(TONE).slice(iR), iR);
         }

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void SphericalHeatAdvection<TONE,TTWO,TTHREE>::sub(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      if(c != 1.0)
      {
         rS.subData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            rS.addSlice((c*r(iR_))*u.comp(TONE).slice(iR), iR);
         }

         rS.subData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.subData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      } else
      {
         rS.subData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            rS.addSlice(r(iR_)*u.comp(TONE).slice(iR), iR);
         }

         rS.subData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

         rS.subData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
      }
   }
}
}

#endif // QUICC_PHYSICAL_SPHERICALHEATADVECTION_HPP
