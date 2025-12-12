/**
 * @file Dot.hpp
 * @brief Implementation of a generic scalar product
 */

#ifndef QUICC_PHYSICAL_DOT_HPP
#define QUICC_PHYSICAL_DOT_HPP

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
    * @brief Implementation of a generic scalar product
    */
   class Dot
   {
      public:
         /**
          * @brief Set S to scalar product
          */
         template <typename TFIELD>
            static void set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add scalar product to S
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract scalar product from S
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

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

      private:
         /// @tparam T scalar
         template <class T = double> struct AddDotFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            AddDotFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            AddDotFunctor() = delete;

            /// @brief dtor
            ~AddDotFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param uk
            /// @param vi
            /// @param vj
            /// @param vk
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T w, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return w + _scaling * (ui * vi + uj * vj + uk * vk);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SubDotFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SubDotFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SubDotFunctor() = delete;

            /// @brief dtor
            ~SubDotFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param uk
            /// @param vi
            /// @param vj
            /// @param vk
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T w, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return w - _scaling * (ui * vi + uj * vj + uk * vk);
            }
         };

         template <typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <typename TFIELD> void Dot::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f)
   {
      for(auto&& [k, flag]: f.enabled())
      {
         vs.push_back(f.comp(k).dataView());
      }
   }

   template <typename TFIELD> void Dot::set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = Pointwise::DotFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, v);
         collectViews(vs, w);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
      {
         auto vIt = v.data().cbegin();
         auto wIt = w.data().cbegin();

         if(c != 1.0)
         {
            rS.setData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
            ++vIt;
            ++wIt;

            for(; vIt != v.data().cend(); ++vIt,++wIt)
            {
               rS.addData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
            }
         } else
         {
            rS.setData((vIt->second.data().array()*wIt->second.data().array()).matrix());
            ++vIt;
            ++wIt;

            for(; vIt != v.data().cend(); ++vIt,++wIt)
            {
               rS.addData((vIt->second.data().array()*wIt->second.data().array()).matrix());
            }
         }
      }
   }

   template <typename TFIELD> void Dot::add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = AddDotFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, v);
         collectViews(vs, w);
         op.apply(rS.rDataView(), rS.dataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
      {
         auto wIt = w.data().cbegin();
         if(c != 1.0)
         {
            for(auto vIt = v.data().cbegin(); vIt != v.data().cend(); ++vIt)
            {
               rS.addData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
               ++wIt;
            }
         } else
         {
            for(auto vIt = v.data().cbegin(); vIt != v.data().cend(); ++vIt)
            {
               rS.addData((vIt->second.data().array()*wIt->second.data().array()).matrix());
               ++wIt;
            }
         }
      }
   }

   template <typename TFIELD> void Dot::sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = SubDotFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vViews;
         std::vector<view_t> vs;
         collectViews(vs, v);
         collectViews(vs, w);
         op.apply(rS.rDataView(), rS.dataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
      {
         auto wIt = w.data().cbegin();
         if(c != 1.0)
         {
            for(auto vIt = v.data().cbegin(); vIt != v.data().cend(); ++vIt)
            {
               rS.subData(c*(vIt->second.data().array()*wIt->second.data().array()).matrix());
               ++wIt;
            }
         } else
         {
            for(auto vIt = v.data().cbegin(); vIt != v.data().cend(); ++vIt)
            {
               rS.subData((vIt->second.data().array()*wIt->second.data().array()).matrix());
               ++wIt;
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_DOT_HPP
