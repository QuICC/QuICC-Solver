/**
 * @file SphericalZComponent.hpp
 * @brief Implementation of the spherical Z component of a field
 */

#ifndef QUICC_PHYSICAL_SPHERICALZCOMPONENT_HPP
#define QUICC_PHYSICAL_SPHERICALZCOMPONENT_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "ViewOps/Slicewise/NoGridOp.hpp"
#include "QuICC/PhysicalOperators/details/FunctorHelpers.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical Z component of a field
    */
   class SphericalZComponent
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD>
            static void set(TFIELD &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void set(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void add(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void sub(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c = 1.0);


      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalZComponent() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalZComponent() = default;

      private:
         /// Functor to map resolution object
         struct IdxResFunctor
         {
            const Resolution& _res;

            IdxResFunctor(const Resolution& res) : _res(res) {};

            /// @brief deleted default constructor
            IdxResFunctor() = delete;

            /// @brief dtor
            ~IdxResFunctor() = default;

            int dim2D(const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(k);
            }

            int idx2D(const int j, const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(j, k);
            }

            int dim3D() const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
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
            /// @param gc
            /// @param gs
            /// @param ui
            /// @param uj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gc, T gs, T ui, T uj)
            {
               return _scaling * (gc * ui - gs * uj);
            }
         };
   };

   template <typename TFIELD>
   void SphericalZComponent::set(TFIELD &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set(rS, f, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD>
   void SphericalZComponent::add(TFIELD &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add(rS, f, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD>
   void SphericalZComponent::sub(TFIELD &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub(rS, f, cosTheta, sinTheta, v, c);
   }



   template <typename TFIELD, typename TIDXFUNC>
   void SphericalZComponent::set(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = SetFunctor<scalar_t>;
         fct_t f(c);
         grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
         grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
         Slicewise::Cpu::NoGridOp<1, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t> op(f);
         op.apply(rS.rGlobalView(), vCos, vSin, v.comp(FieldComponents::Physical::R).globalView(), v.comp(FieldComponents::Physical::THETA).globalView());
      }
      else
      {
         int nR = idxFunc.dim3D();
         int nTh;
         int iTh_;

         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.setProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  rS.subProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               }
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.setProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  rS.subProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalZComponent::add(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = details::AddTmplFunctor<scalar_t, SetFunctor>;
         fct_t f(c);
         grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
         grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
         Slicewise::Cpu::NoGridOp<1, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rGlobalView(), vCos, vSin, v.comp(FieldComponents::Physical::R).globalView(), v.comp(FieldComponents::Physical::THETA).globalView(), rS.dataView());
      }
      else
      {
         int nR = idxFunc.dim3D();
         int nTh;
         int iTh_;

         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.addProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  rS.subProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               }
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.addProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  rS.subProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalZComponent::sub(TFIELD &rS, const TIDXFUNC& idxFunc, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = details::SubTmplFunctor<scalar_t, SetFunctor>;
         fct_t f(c);
         grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
         grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
         Slicewise::Cpu::NoGridOp<1, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t> op(f);
         op.apply(rS.rGlobalView(), vCos, vSin, v.comp(FieldComponents::Physical::R).globalView(), v.comp(FieldComponents::Physical::THETA).globalView(), rS.dataView());
      }
      else
      {
         int nR = idxFunc.dim3D();
         int nTh;
         int iTh_;

         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.subProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
                  rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               }
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.subProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
                  rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               }
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALZCOMPONENT_HPP
