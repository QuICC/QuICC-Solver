/**
 * @file Dot.hpp
 * @brief Implementation of a generic scalar product
 */

#ifndef QUICC_PHYSICAL_DOT_HPP
#define QUICC_PHYSICAL_DOT_HPP

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
    * @brief Implementation of a generic scalar product
    */
   class Dot
   {
      public:
         /**
          * @brief Set S to scalar product
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add scalar product to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract scalar product from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         Dot() = default;

         /**
          * @brief Empty destructor
          */
         ~Dot() = default;
   };

   inline void Dot::set(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      auto vIt = v.data().cbegin();
      auto wIt = w.data().cbegin();

      if(c != 1.0)
      {
         rS.setData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
         ++vIt;
         ++wIt;

         for(; vIt != v.data.cend(); ++vIt,++wIt)
         {
            rS.addData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
         }
      } else
      {
         rS.setData((vIt->second.data().array()*wIt->second.data().array()).matrix());
         ++vIt;
         ++wIt;

         for(; vIt != v.data.cend(); ++vIt,++wIt)
         {
            rS.addData((vIt->second.data().array()*wIt->second.data().array()).matrix());
         }
      }
   }

   inline void Dot::add(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      auto wIt = w.data().cbegin();
      if(c != 1.0)
      {
         for(auto vIt = v.data().cbegin(); vIt != v.data.cend(); ++vIt)
         {
            rS.addData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
            ++wIt;
         }
      } else
      {
         for(auto vIt = v.data().cbegin(); vIt != v.data.cend(); ++vIt)
         {
            rS.addData((vIt->second.data().array()*wIt->second.data().array()).matrix());
            ++wIt;
         }
      }
   }

   inline void Dot::sub(Framework::Selector::PhysicalScalarField &rS, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      auto wIt = w.data().cbegin();
      if(c != 1.0)
      {
         for(auto vIt = v.data().cbegin(); vIt != v.data.cend(); ++vIt)
         {
            rS.subData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
            ++wIt;
         }
      } else
      {
         for(auto vIt = v.data().cbegin(); vIt != v.data.cend(); ++vIt)
         {
            rS.subData((vIt->second.data().array()*wIt->second.data().array()).matrix());
            ++wIt;
         }
      }
   }
}
}

#endif // QUICC_PHYSICAL_DOT_HPP
