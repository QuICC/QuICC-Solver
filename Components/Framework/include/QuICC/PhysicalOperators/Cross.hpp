/**
 * @file Cross.hpp
 * @brief Implementation of a generic vector cross product
 */

#ifndef QUICC_PHYSICAL_CROSS_HPP
#define QUICC_PHYSICAL_CROSS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "ViewOps/Pointwise/Functors.hpp"
#include "ViewOps/Pointwise/Pointwise.hpp"

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
         template <typename TFIELD>
            static void set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add cross product component to S
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract cross product component from S
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

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

      private:
         /// @tparam T scalar
         template <class T = double> struct AddCrossCompFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            AddCrossCompFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            AddCrossCompFunctor() = delete;

            /// @brief dtor
            ~AddCrossCompFunctor() = default;

            /// @brief Cross product, component wise
            /// @param uj
            /// @param uk
            /// @param vj
            /// @param vk
            /// @return i component of cross product
            QUICC_CUDA_HOSTDEV T operator()(T w, T uj, T uk, T vj, T vk)
            {
               return w + _scaling * (uj * vk - uk * vj);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SubCrossCompFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SubCrossCompFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SubCrossCompFunctor() = delete;

            /// @brief dtor
            ~SubCrossCompFunctor() = default;

            /// @brief Cross product, component wise
            /// @param uj
            /// @param uk
            /// @param vj
            /// @param vk
            /// @return i component of cross product
            QUICC_CUDA_HOSTDEV T operator()(T w, T uj, T uk, T vj, T vk)
            {
               return w - _scaling * (uj * vk - uk * vj);
            }
         };
   };

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> template <typename TFIELD> inline void Cross<TFIRST,TSECOND>::set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = Pointwise::CrossCompFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rDataView(), v.comp(TFIRST).dataView(), v.comp(TSECOND).dataView(), w.comp(TFIRST).dataView(), w.comp(TSECOND).dataView());
      }
      else
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
   }

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> template <typename TFIELD> inline void Cross<TFIRST,TSECOND>::add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = AddCrossCompFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rDataView(), rS.dataView(), v.comp(TFIRST).dataView(), v.comp(TSECOND).dataView(), w.comp(TFIRST).dataView(), w.comp(TSECOND).dataView());
      }
      else
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
   }

   template <FieldComponents::Physical::Id TFIRST,FieldComponents::Physical::Id TSECOND> template <typename TFIELD> inline void Cross<TFIRST,TSECOND>::sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = SubCrossCompFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rDataView(), rS.dataView(), v.comp(TFIRST).dataView(), v.comp(TSECOND).dataView(), w.comp(TFIRST).dataView(), w.comp(TSECOND).dataView());
      }
      else
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
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_CROSS_HPP
