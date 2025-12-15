/**
 * @file SphericalBuoyancy.hpp
 * @brief Implementation of the spherical buoyancy term
 */

#ifndef QUICC_PHYSICAL_SPHERICALBUOYANCY_HPP
#define QUICC_PHYSICAL_SPHERICALBUOYANCY_HPP

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
    * @brief Implementation of the spherical buoyancy term
    */
   class SphericalBuoyancy
   {
      public:
         /**
          * @brief Set S
          */
         template <typename TFIELD>
            static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const TFIELD &q, const MHDFloat c = 1.0);

         /**
          * @brief Add to S
          */
         template <typename TFIELD>
            static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const TFIELD &q, const MHDFloat c = 1.0);

         /**
          * @brief Substract from S
          */
         template <typename TFIELD>
            static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const TFIELD &v, const MHDFloat c = 1.0);

         /**
          * @brief Set S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& r, const TFIELD &q, const MHDFloat c = 1.0);

         /**
          * @brief Add to S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& r, const TFIELD &q, const MHDFloat c = 1.0);

         /**
          * @brief Substract from S
          */
         template <typename TFIELD, typename TIDXFUNC>
            static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& r, const TFIELD &v, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalBuoyancy() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalBuoyancy() = default;

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

            int dim3D() const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
            }

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
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T g, T ui)
            {
               return _scaling * (g * ui);
            }
         };
   };

   template <typename TFIELD>
      void SphericalBuoyancy::set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const TFIELD &q, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set(rS, compId, f, r, q, c);
   }

   template <typename TFIELD>
      void SphericalBuoyancy::add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const TFIELD &q, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add(rS, compId, f, r, q, c);
   }

   template <typename TFIELD>
      void SphericalBuoyancy::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const TFIELD &q, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub(rS, compId, f, r, q, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
      void SphericalBuoyancy::set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& r, const TFIELD &q, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetFunctor<scalar_t>;
            fct_t f(c);
            grid_t vGrid(const_cast<scalar_t *>(r.data()), r.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vGrid, q.globalView());
         }
         else
         {
            int nR = idxFunc.dim3D();
            int iR_;
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  rS.setSlice(c*r(iR_)*q.slice(iR), iR);
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  rS.setSlice(r(iR_)*q.slice(iR), iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         rS.setZeros();
      } else if(compId == FieldComponents::Physical::PHI)
      {
         rS.setZeros();
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
      void SphericalBuoyancy::add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& r, const TFIELD &q, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetFunctor>;
            fct_t f(c);
            grid_t vGrid(const_cast<scalar_t *>(r.data()), r.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vGrid, q.globalView(), rS.dataView());
         }
         else
         {
            int nR = idxFunc.dim3D();
            int iR_;

            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  rS.addSlice(c*r(iR_)*q.slice(iR), iR);
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  rS.addSlice(r(iR_)*q.slice(iR), iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
      void SphericalBuoyancy::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& r, const TFIELD &q, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetFunctor>;
            fct_t f(c);
            grid_t vGrid(const_cast<scalar_t *>(r.data()), r.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vGrid, q.globalView(), rS.dataView());
         }
         else
         {
            int nR = idxFunc.dim3D();
            int iR_;
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  rS.subSlice(c*r(iR_)*q.slice(iR), iR);
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  rS.subSlice(r(iR_)*q.slice(iR), iR);
               }
            }
         }
      }
   }

} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALBUOYANCY_HPP
