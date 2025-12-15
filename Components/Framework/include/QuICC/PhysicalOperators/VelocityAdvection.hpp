/**
 * @file VelocityAdvection.hpp
 * @brief Implementation of a generic primitive velocity advection
 */

#ifndef QUICC_PHYSICAL_VELOCITYADVECTION_HPP
#define QUICC_PHYSICAL_VELOCITYADVECTION_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "ViewOps/Pointwise/Pointwise.hpp"
#include "QuICC/PhysicalOperators/details/FunctorHelpers.hpp"

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
         template <typename TFIELD>
            static void set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

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

      private:
         /// @tparam T scalar
         template <class T = double> struct VelAdvFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            VelAdvFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            VelAdvFunctor() = delete;

            /// @brief dtor
            ~VelAdvFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param uk
            /// @param vi
            /// @param vj
            /// @param vk
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return _scaling * (ui * vi + uj * vj + uk * vk);
            }
         };

         template <typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,TTHREE>::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD &f)
      {
         vs.push_back(f.comp(TONE).dataView());
         vs.push_back(f.comp(TTWO).dataView());
         vs.push_back(f.comp(TTHREE).dataView());
      }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,TTHREE>::set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = VelAdvFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
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
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,TTHREE>::add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = details::AddTmplFunctor<scalar_t, VelAdvFunctor>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5), rS.dataView());
      }
      else
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
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,TTHREE>::sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = details::SubTmplFunctor<scalar_t, VelAdvFunctor>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5), rS.dataView());
      }
      else
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
         template <typename TFIELD>
            static void set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add primitive velocity advection product to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract primitive velocity advection product from S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)q = u_x\partial_x q + u_y\partial_y q + u_z\partial_z q\f$
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

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

      private:
         /// @tparam T scalar
         template <class T = double> struct VelAdvFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            VelAdvFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            VelAdvFunctor() = delete;

            /// @brief dtor
            ~VelAdvFunctor() = default;

            /// @brief Dot product
            /// @param ui
            /// @param uj
            /// @param vi
            /// @param vj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T ui, T uj, T vi, T vj)
            {
               return _scaling * (ui * vi + uj * vj);
            }
         };

         template <typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD &f)
      {
         vs.push_back(f.comp(TONE).dataView());
         vs.push_back(f.comp(TTWO).dataView());
      }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::set(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = VelAdvFunctor<scalar_t>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3));
      }
      else
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
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::add(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = details::AddTmplFunctor<scalar_t, VelAdvFunctor>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), rS.dataView());
      }
      else
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
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO>
      template <typename TFIELD>
      void VelocityAdvection<TONE,TTWO,FieldComponents::Physical::NOTUSED>::sub(TFIELD &rS, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using fct_t = details::SubTmplFunctor<scalar_t, VelAdvFunctor>;
         fct_t f(c);
         Pointwise::Cpu::Op<fct_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         op.apply(rS.rDataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), rS.dataView());
      }
      else
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
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_VELOCITYADVECTION_HPP
