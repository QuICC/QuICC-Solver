/**
 * @file StreamAdvection.hpp
 * @brief Implementation of a generic streamfunction advection
 */

#ifndef QUICC_PHYSICAL_STREAMADVECTION_HPP
#define QUICC_PHYSICAL_STREAMADVECTION_HPP

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
    * @brief Implementation of a generic streamfunction advection
    */
   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> class StreamAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)q = -\partial_y\psi\partial_x q + \partial_x\psi\partial_y q\f$
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)q = -\partial_y\psi\partial_x q + \partial_x\psi\partial_y q\f$
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)q = -\partial_y\psi\partial_x q + \partial_x\psi\partial_y q\f$
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         StreamAdvection() = default;

         /**
          * @brief Empty destructor
          */
         ~StreamAdvection() = default;
   };

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamAdvection<TXComp,TYComp>::set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array()).matrix());
      } else
      {
         rS.setData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamAdvection<TXComp,TYComp>::add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamAdvection<TXComp,TYComp>::sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.addData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.addData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array()).matrix());
      }
   }
}
}

#endif // QUICC_PHYSICAL_STREAMADVECTION_HPP
