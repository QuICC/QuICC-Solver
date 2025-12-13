/**
 * @file SphericalHeatAdvection.hpp
 * @brief Implementation of a spherical internal heat advection
 */

#ifndef QUICC_PHYSICAL_SPHERICALHEATADVECTION_HPP
#define QUICC_PHYSICAL_SPHERICALHEATADVECTION_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "ViewOps/Slicewise/NoGridOp.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of a spherical internal heating advection
    */
   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE> class SphericalHeatAdvection
   {
      public:
         /**
          * @brief Set S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q+q_b)\f$
          */
         template <typename TFIELD>
          static void set(TFIELD &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q + q_b)\f$
          */
         template <typename TFIELD>
          static void add(TFIELD &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q + q_b)\f$
          */
         template <typename TFIELD>
          static void sub(TFIELD &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Set S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q+q_b)\f$
          */
         template <typename TFIELD, typename TIDXFUNC>
          static void set(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Add to S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q + q_b)\f$
          */
         template <typename TFIELD, typename TIDXFUNC>
          static void add(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

         /**
          * @brief Substract S
          *
          *    \f$ \left(\vec u\cdot\nabla\right)(q + q_b)\f$
          */
         template <typename TFIELD, typename TIDXFUNC>
          static void sub(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalHeatAdvection() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalHeatAdvection() = default;

      private:
         /// Functor to map resolution object
         struct IdxResFunctor
         {
            const Resolution& _res;

            const int dim3D;

            IdxResFunctor(const Resolution& res) : _res(res), dim3D(res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>()) {};

            /// @brief deleted default constructor
            IdxResFunctor() = delete;

            /// @brief dtor
            ~IdxResFunctor() = default;

            int idx3D(const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetFunctor() = delete;

            /// @brief dtor
            ~SetFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @param uj
            /// @param uk
            /// @param vi
            /// @param vj
            /// @param vk
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return _scaling * (ui * vi + uj * vj + uk * vk - g * ui);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct AddFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            AddFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            AddFunctor() = delete;

            /// @brief dtor
            ~AddFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @param uj
            /// @param uk
            /// @param vi
            /// @param vj
            /// @param vk
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T w, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return w + _scaling * (ui * vi + uj * vj + uk * vk - g * ui);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SubFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SubFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SubFunctor() = delete;

            /// @brief dtor
            ~SubFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @param uj
            /// @param uk
            /// @param vi
            /// @param vj
            /// @param vk
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T w, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return w - _scaling * (ui * vi + uj * vj + uk * vk - g * ui);
            }
         };

         template <typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
      template <typename TFIELD>
      void SphericalHeatAdvection<TONE,TTWO,TTHREE>::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD &f)
      {
         vs.push_back(f.comp(TONE).dataView());
         vs.push_back(f.comp(TTWO).dataView());
         vs.push_back(f.comp(TTHREE).dataView());
      }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
   template <typename TFIELD>
   void SphericalHeatAdvection<TONE,TTWO,TTHREE>::set(TFIELD &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set<TONE,TTWO,TTHREE>(rS, f, r, u, gradQ, c);
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
   template <typename TFIELD>
   void SphericalHeatAdvection<TONE,TTWO,TTHREE>::add(TFIELD &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add<TONE,TTWO,TTHREE>(rS, f, r, u, gradQ, c);
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
   template <typename TFIELD>
   void SphericalHeatAdvection<TONE,TTWO,TTHREE>::sub(TFIELD &rS, const Resolution& res, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub<TONE,TTWO,TTHREE>(rS, f, r, u, gradQ, c);
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
   template <typename TFIELD, typename TIDXFUNC>
   void SphericalHeatAdvection<TONE,TTWO,TTHREE>::set(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = SetFunctor<scalar_t>;
         fct_t f(c);
         grid_t vGrid(const_cast<scalar_t *>(r.data()), r.size());
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         Slicewise::Cpu::NoGridOp<2, fct_t, view_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rGlobalView(), vGrid, vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
      {
         int nR = idxFunc.dim3D;
         int iR_;

         if(c != 1.0)
         {
            rS.setData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               rS.subSlice((c*r(iR_))*u.comp(TONE).slice(iR), iR);
            }

            rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

            rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
         } else
         {
            rS.setData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               rS.subSlice(r(iR_)*u.comp(TONE).slice(iR), iR);
            }

            rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

            rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
         }
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
   template <typename TFIELD, typename TIDXFUNC>
   void SphericalHeatAdvection<TONE,TTWO,TTHREE>::add(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = AddFunctor<scalar_t>;
         fct_t f(c);
         grid_t vGrid(const_cast<scalar_t *>(r.data()), r.size());
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         Slicewise::Cpu::NoGridOp<2, fct_t, view_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rGlobalView(), vGrid, rS.dataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
      {
         int nR = idxFunc.dim3D;
         int iR_;

         if(c != 1.0)
         {
            rS.addData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               rS.subSlice((c*r(iR_))*u.comp(TONE).slice(iR), iR);
            }

            rS.addData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

            rS.addData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
         } else
         {
            rS.addData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               rS.subSlice(r(iR_)*u.comp(TONE).slice(iR), iR);
            }

            rS.addData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

            rS.addData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
         }
      }
   }

   template <FieldComponents::Physical::Id TONE, FieldComponents::Physical::Id TTWO, FieldComponents::Physical::Id TTHREE>
   template <typename TFIELD, typename TIDXFUNC>
   void SphericalHeatAdvection<TONE,TTWO,TTHREE>::sub(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& r, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &u, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &gradQ, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = SubFunctor<scalar_t>;
         fct_t f(c);
         grid_t vGrid(const_cast<scalar_t *>(r.data()), r.size());
         std::vector<view_t> vs;
         collectViews(vs, u);
         collectViews(vs, gradQ);
         Slicewise::Cpu::NoGridOp<2, fct_t, view_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rGlobalView(), vGrid, rS.dataView(), vs.at(0), vs.at(1), vs.at(2), vs.at(3), vs.at(4), vs.at(5));
      }
      else
      {
         int nR = idxFunc.dim3D;
         int iR_;

         if(c != 1.0)
         {
            rS.subData(c*(u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               rS.addSlice((c*r(iR_))*u.comp(TONE).slice(iR), iR);
            }

            rS.subData(c*(u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

            rS.subData(c*(u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
         } else
         {
            rS.subData((u.comp(TONE).data().array()*gradQ.comp(TONE).data().array()).matrix());
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               rS.addSlice(r(iR_)*u.comp(TONE).slice(iR), iR);
            }

            rS.subData((u.comp(TTWO).data().array()*gradQ.comp(TTWO).data().array()).matrix());

            rS.subData((u.comp(TTHREE).data().array()*gradQ.comp(TTHREE).data().array()).matrix());
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALHEATADVECTION_HPP
