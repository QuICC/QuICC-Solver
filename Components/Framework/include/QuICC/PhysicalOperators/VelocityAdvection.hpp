/**
 * @file VelocityAdvection.hpp
 * @brief Implementation of a generic primitive velocity advection
 */

#ifndef QUICC_PHYSICAL_VELOCITYADVECTION_HPP
#define QUICC_PHYSICAL_VELOCITYADVECTION_HPP

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
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of a generic primitive velocity advection in 3D space
    */
   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> class VelocityAdvection
   {
      public:
         /**
          * @brief Set S to primitive velocity advection product
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         VelocityAdvection() = default;

         /**
          * @brief Empty destructor
          */
         ~VelocityAdvection() = default;
   };

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void VelocityAdvection<TONE,TTWO,TTHREE>::set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
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

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void VelocityAdvection<TONE,TTWO,TTHREE>::add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
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

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> void VelocityAdvection<TONE,TTWO,TTHREE>::sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
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

   /**
    * @brief Implementation of a generic primitive velocity advection in 2D space
    */
   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> class VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>
   {
      public:
         /**
          * @brief Set S to primitive velocity advection product
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
          static void sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         VelocityAdvection() = default;

         /**
          * @brief Empty destructor
          */
         ~VelocityAdvection() = default;
   };

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());
      } else
      {
         rS.setData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());
      } else
      {
         rS.addData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO> void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.subData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());
      } else
      {
         rS.subData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());

         rS.subData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());
      }
   }
}
}

#endif // QUICC_PHYSICAL_VELOCITYADVECTION_HPP
