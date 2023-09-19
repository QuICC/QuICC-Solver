/**
 * @file Cross.hpp
 * @brief Implementation of a generic vector cross product
 */

#ifndef QUICC_PHYSICAL_CROSS_HPP
#define QUICC_PHYSICAL_CROSS_HPP

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
    * @brief Implementation of a generic vector cross product
    */
   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> class Cross
   {
      public:
         /**
          * @brief Set S to cross product component
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add cross product component to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract cross product component from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         Cross() = default;

         /**
          * @brief Empty destructor
          */
         ~Cross() = default;
   };

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> inline void Cross<TFIRST,TSECOND>::set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.setData(c*(v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData(c*(v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      } else
      {
         rS.setData((v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData((v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> inline void Cross<TFIRST,TSECOND>::add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.addData(c*(v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData(c*(v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      } else
      {
         rS.addData((v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.subData((v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      }
   }

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> inline void Cross<TFIRST,TSECOND>::sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      if(c != 1.0)
      {
         rS.subData(c*(v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.addData(c*(v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      } else
      {
         rS.subData((v.comp(TFIRST).data().array() * w.comp(TSECOND).data().array()).matrix());

         rS.addData((v.comp(TSECOND).data().array() * w.comp(TFIRST).data().array()).matrix());
      }
   }
}
}

#endif // QUICC_PHYSICAL_CROSS_HPP
