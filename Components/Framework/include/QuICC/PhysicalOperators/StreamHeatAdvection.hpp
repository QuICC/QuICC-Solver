/**
 * @file StreamHeatAdvection.hpp
 * @brief Implementation of a generic streamfunction advection including conducting state
 */

#ifndef QUICC_PHYSICAL_STREAMHEATADVECTION_HPP
#define QUICC_PHYSICAL_STREAMHEATADVECTION_HPP

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
    * @brief Implementation of a generic streamfunction advection including conducting state
    */
   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> class StreamHeatAdvection
   {
      public:
         /**
          * @brief Set S to streamfunction advection product
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         StreamHeatAdvection() = default;

         /**
          * @brief Empty destructor
          */
         ~StreamHeatAdvection() = default;
   };

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamHeatAdvection<TXComp,TYComp>::set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      } else
      {
         rS.setData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamHeatAdvection<TXComp,TYComp>::add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      } else
      {
         rS.addData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.subData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp> void StreamHeatAdvection<TXComp,TYComp>::sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.addData(c*(dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      } else
      {
         rS.subData((dPsi.comp(TXComp).data().array()*w.comp(TYComp).data().array()).matrix());

         rS.addData((dPsi.comp(TYComp).data().array()*w.comp(TXComp).data().array() + dPsi.comp(TYComp).data().array()).matrix());
      }
   }
}
}

#endif // QUICC_PHYSICAL_STREAMHEATADVECTION_HPP
