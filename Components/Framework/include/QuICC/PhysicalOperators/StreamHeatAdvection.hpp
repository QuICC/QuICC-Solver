/**
 * @file StreamHeatAdvection.hpp
 * @brief Implementation of a generic streamfunction advection including conducting state
 */

#ifndef QUICC_PHYSICAL_STREAMHEATADVECTION_HPP
#define QUICC_PHYSICAL_STREAMHEATADVECTION_HPP

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
         template <typename TFIELD>
            static void set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Add streamfunction advection product to S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

         /**
          * @brief Substract streamfunction advection product from S
          *
          *    \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\overline{T} = -\partial_y\psi\partial_x \theta -\partial_y\psi\partial_x x + \partial_x\psi\partial_y \theta\f$
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c = 1.0);

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

      private:
         /// @tparam T scalar
         template <class T = double> struct StreamHeatAdvFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            StreamHeatAdvFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            StreamHeatAdvFunctor() = delete;

            /// @brief dtor
            ~StreamHeatAdvFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param vi
            /// @param vj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T ui, T uj, T vi, T vj)
            {
               return _scaling * (ui * vj - uj * vi - uj);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct AddStreamHeatAdvFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            AddStreamHeatAdvFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            AddStreamHeatAdvFunctor() = delete;

            /// @brief dtor
            ~AddStreamHeatAdvFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param vi
            /// @param vj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T w, T ui, T uj, T vi, T vj)
            {
               return w + _scaling * (ui * vj - uj * vi - uj);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SubStreamHeatAdvFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SubStreamHeatAdvFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SubStreamHeatAdvFunctor() = delete;

            /// @brief dtor
            ~SubStreamHeatAdvFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param vi
            /// @param vj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T w, T ui, T uj, T vi, T vj)
            {
               return w - _scaling * (ui * vj - uj * vi - uj);
            }
         };

         template <typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp>
      template <typename TFIELD>
      void StreamHeatAdvection<TXComp,TYComp>::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD &f)
      {
         vs.push_back(f.comp(TXComp).dataView());
         vs.push_back(f.comp(TYComp).dataView());
      }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp>
      template <typename TFIELD>
      void StreamHeatAdvection<TXComp,TYComp>::set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = StreamHeatAdvFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, dPsi);
         collectViews(vs, w);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3));
      }
      else
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
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp>
      template <typename TFIELD>
      void StreamHeatAdvection<TXComp,TYComp>::add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = AddStreamHeatAdvFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, dPsi);
         collectViews(vs, w);
         op.apply(rS.rDataView(), rS.dataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3));
      }
      else
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
   }

   template <FieldComponents::Physical::Id TXComp, FieldComponents::Physical::Id TYComp>
      template <typename TFIELD>
      void StreamHeatAdvection<TXComp,TYComp>::sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &dPsi, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = SubStreamHeatAdvFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, dPsi);
         collectViews(vs, w);
         op.apply(rS.rDataView(), rS.dataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3));
      }
      else
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
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_STREAMHEATADVECTION_HPP
