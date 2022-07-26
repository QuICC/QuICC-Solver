/**
 * @file FieldTools.hpp
 * @brief Some useful tools to work with fields
 */

#ifndef QUICC_DATATYPES_FIELDTOOLS_HPP
#define QUICC_DATATYPES_FIELDTOOLS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>
#include <tuple>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Single configuration class for the different scalar fields
    */
   class FieldTools
   {
      public:
         /**
          * @brief Create const field data information
          */
         template <typename T, template <typename> class TField> static std::vector<std::tuple<int, int , const T *> > createInfo(const TField<T>& field);

         /**
          * @brief Create field data information
          */
         template <typename T, template <typename> class TField> static std::vector<std::tuple<int, int , T *> >  createInfo(TField<T>& rField);

         /**
          * @brief Combine two field with arithmetic operation
          */
         template <typename T1, typename T2, template <typename> class TField> static void combine(TField<T1>& rRhs, const TField<T2>& lhs, const std::size_t arithId);

         /**
          * @brief Set field to negative
          */
         template <typename T, template <typename> class TField> static void negative(TField<T>& rRhs);

      protected:

      private:
         /**
          * @brief Constructor
          */
         FieldTools();

         /**
          * @brief Destructor
          */
         virtual ~FieldTools();
   };

   template <typename T, template <typename> class TField> std::vector<std::tuple<int, int , const T *> > FieldTools::createInfo(const TField<T>& field)
   {
      // Storage for the field information
      std::vector<std::tuple<int,int, const T *> > fieldInfo;

      // Loop over all slices
      fieldInfo.reserve(field.nSlice());
      for(int i = 0; i < field.nSlice(); i++)
      {
         int rows = field.slice(i).rows();
         int cols = field.slice(i).cols();

         fieldInfo.push_back(std::make_tuple(rows, cols, field.data(i)));
      }

      return fieldInfo;
   }

   template <typename T, template <typename> class TField> std::vector<std::tuple<int, int , T *> > FieldTools::createInfo(TField<T>& rField)
   {
      // Storage for the field information
      std::vector<std::tuple<int,int, T *> > fieldInfo;

      // Loop over all slices
      fieldInfo.reserve(rField.nSlice());
      for(int i = 0; i < rField.nSlice(); i++)
      {
         int rows = rField.slice(i).rows();
         int cols = rField.slice(i).cols();

         fieldInfo.push_back(std::make_tuple(rows, cols, rField.rData(i)));
      }

      return fieldInfo;
   }

   template <typename T1, typename T2, template <typename> class TField> void FieldTools::combine(TField<T1>& rRhs, const TField<T2>& lhs, const std::size_t arithId)
   {
      Profiler::RegionFixture<3> fix("FieldTools::combine");
      assert(arithId == Arithmetics::Add::id() || arithId == Arithmetics::Sub::id());

      if constexpr(std::is_same<T1,T2>::value)
      {
         if(arithId == Arithmetics::Add::id())
         {
            rRhs.addData(lhs.data());
         } else if(arithId == Arithmetics::Sub::id())
         {
            rRhs.subData(lhs.data());
         } else
         {
            throw std::logic_error("Unknown arithmetic operation in combine!");
         }
      } else
      {
         throw std::logic_error("Types are incompatible for arithmetic operation in combine");
      }
   }

   template <typename T, template <typename> class TField> void FieldTools::negative(TField<T>& rRhs)
   {
      Profiler::RegionFixture<3> fix("FieldTools::negative");
      rRhs.rData() = -rRhs.data();
   }

}
}

#endif // QUICC_DATATYPES_FIELDTOOLS_HPP
