/*
 * @file SphericalPoincare.hpp
 * @brief Implementation of the spherical poincare term
 */

#ifndef QUICC_PHYSICAL_SPHERICALPOINCARE_HPP
#define QUICC_PHYSICAL_SPHERICALPOINCARE_HPP

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
    * @brief Implementation of the spherical poincare term
    */
   class SphericalPoincare
   {
      public:
         /**
          * @brief Set S to Poincare term
          */
         template <typename TFIELD>
         static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Add Poincare term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Substract Poincare term from S
          */
         template <typename TFIELD>
         static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Set S to Poincare term
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Add Poincare term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Substract Poincare term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPoincare() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalPoincare() = default;

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

            int idx3D(const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetTFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetTFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetTFunctor() = delete;

            /// @brief dtor
            ~SetTFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gR, T gP)
            {
               return -_scaling * gR * gP;
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetPFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetPFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetPFunctor() = delete;

            /// @brief dtor
            ~SetPFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gR, T gT, T gP)
            {
               return _scaling * gR * gT * gP;
            }
         };

   };

   template <typename TFIELD>
   void SphericalPoincare::set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set(rS, compId, f, rGrid, thGrid, phGrid, t, alpha, c);
   }

   template <typename TFIELD>
   void SphericalPoincare::add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add(rS, compId, f, rGrid, thGrid, phGrid, t, alpha, c);
   }

   template <typename TFIELD>
   void SphericalPoincare::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub(rS, compId, f, rGrid, thGrid, phGrid, t, alpha, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPoincare::set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         rS.setZeros();
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetTFunctor<scalar_t>;
            fct_t f(csa);
            Array cos_pt = (phGrid.array() + t).array().cos();
            grid_t vR(const_cast<scalar_t *>(rGrid.data()), rGrid.size());
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            Slicewise::Cpu::NoGridOp<5, fct_t, view_t, 1, 1, 0, grid_t, grid_t> op(f);
            op.apply(rS.rGlobalView(), vR, vCosPt);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               coeff = -csa*rGrid(idxFunc.idx3D(iR));
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  rS.setProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetPFunctor<scalar_t>;
            fct_t f(csa);
            Array cos_T = thGrid.array().array().cos();
            Array sin_Pt = (phGrid.array() + t).array().sin();
            grid_t vR(const_cast<scalar_t *>(rGrid.data()), rGrid.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_Pt.data()), sin_Pt.size());
            Slicewise::Cpu::NoGridOp<10, fct_t, view_t, 1, 1, 1, grid_t, grid_t, grid_t> op(f);
            op.apply(rS.rGlobalView(), vR, vCosT, vSinPt);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               coeff = csa*rGrid(idxFunc.idx3D(iR));
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.setProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPoincare::add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         //
         // Zero
         //
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetTFunctor>;
            fct_t f(csa);
            Array cos_pt = (phGrid.array() + t).array().cos();
            grid_t vR(const_cast<scalar_t *>(rGrid.data()), rGrid.size());
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            Slicewise::Cpu::NoGridOp<5, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vR, vCosPt, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               coeff = csa*rGrid(idxFunc.idx3D(iR));
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  rS.subProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetPFunctor>;
            fct_t f(csa);
            Array cos_T = thGrid.array().array().cos();
            Array sin_Pt = (phGrid.array() + t).array().sin();
            grid_t vR(const_cast<scalar_t *>(rGrid.data()), rGrid.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_Pt.data()), sin_Pt.size());
            Slicewise::Cpu::NoGridOp<10, fct_t, view_t, 1, 1, 1, grid_t, grid_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vR, vCosT, vSinPt, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               coeff = csa*rGrid(idxFunc.idx3D(iR));
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.addProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPoincare::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         //
         // Zero
         //
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetTFunctor>;
            fct_t f(csa);
            Array cos_pt = (phGrid.array() + t).array().cos();
            grid_t vR(const_cast<scalar_t *>(rGrid.data()), rGrid.size());
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            Slicewise::Cpu::NoGridOp<5, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vR, vCosPt, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               coeff = csa*rGrid(idxFunc.idx3D(iR));
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetPFunctor>;
            fct_t f(csa);
            Array cos_T = thGrid.array().array().cos();
            Array sin_Pt = (phGrid.array() + t).array().sin();
            grid_t vR(const_cast<scalar_t *>(rGrid.data()), rGrid.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_Pt.data()), sin_Pt.size());
            Slicewise::Cpu::NoGridOp<10, fct_t, view_t, 1, 1, 1, grid_t, grid_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vR, vCosT, vSinPt, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               coeff = csa*rGrid(idxFunc.idx3D(iR));
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.subProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
               }
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALPOINCARE_HPP
