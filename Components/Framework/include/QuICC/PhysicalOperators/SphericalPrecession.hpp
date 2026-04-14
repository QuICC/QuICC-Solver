/*
 * @file SphericalPrecession.hpp
 * @brief Implementation of the spherical Coriolis + precession term
 */

#ifndef QUICC_PHYSICAL_SPHERICALPRECESSION_HPP
#define QUICC_PHYSICAL_SPHERICALPRECESSION_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/PhysicalOperators/details/FunctorHelpers.hpp"
#include "ViewOps/Slicewise/Cpu/NoGridOp.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical Coriolis + precession term
    */
   class SphericalPrecession
   {
      public:
         /**
          * @brief Set S to Coriolis + precession term
          */
         template <typename TFIELD>
         static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis + precesion term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Substract (Coriolis + precession) term from S
          */
         template <typename TFIELD>
         static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Set S to Coriolis + precession term
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis + precesion term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Substract (Coriolis + precession) term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPrecession() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalPrecession() = default;

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

            int dim3D() const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
            }

            int idx2D(const int j, const int k) const
            {
               return _res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(j, k);
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetRFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _sA;
            T _sB;

            /// @brief ctor
            /// @param scaling
            SetRFunctor(T sA, T sB) : _sA(sA), _sB(sB) {};

            /// @brief deleted default constructor
            SetRFunctor() = delete;

            /// @brief dtor
            ~SetRFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gcT, T gsT, T gcP, T gsP, T ui, T uj)
            {
               return _sA * (gsP * ui + gcT * gcP * uj) - _sB * gsT * uj;
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetTFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _sA;
            T _sB;

            /// @brief ctor
            /// @param scaling
            SetTFunctor(T sA, T sB) : _sA(sA), _sB(sB) {};

            /// @brief deleted default constructor
            SetTFunctor() = delete;

            /// @brief dtor
            ~SetTFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gcT, T gsT, T gcP, T gsP, T ui, T uj)
            {
               return -_sA * (gsP * ui + gsT * gcP * uj) - _sB * gcT * uj;
            }
         };

         /// @tparam T scalar
         template <class T = double> struct SetPFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _sA;
            T _sB;

            /// @brief ctor
            /// @param scaling
            SetPFunctor(T sA, T sB) : _sA(sA), _sB(sB) {};

            /// @brief deleted default constructor
            SetPFunctor() = delete;

            /// @brief dtor
            ~SetPFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gcT, T gsT, T gcP, T gsP, T ui, T uj)
            {
               return _sA * gcP * (-gcT * ui + gsT * uj) + _sB * (gsT * ui + gcT * uj);
            }
         };
   };

   template <typename TFIELD>
   void SphericalPrecession::set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set(rS, compId, f, thGrid, phGrid, v, t, alpha, corC, preC, c);
   }

   template <typename TFIELD>
   void SphericalPrecession::add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add(rS, compId, f, thGrid, phGrid, v, t, alpha, corC, preC, c);
   }

   template <typename TFIELD>
   void SphericalPrecession::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub(rS, compId, f, thGrid, phGrid, v, t, alpha, corC, preC, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPrecession::set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + preC*std::cos(alpha));
      int nR = idxFunc.dim3D();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetRFunctor<scalar_t>;
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::THETA).dataView();
            auto vP = v.comp(FieldComponents::Physical::PHI).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // Theta component
                  rS.setProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
                  // Phi components
                  coeff = cA*std::cos(theta);
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::sin(theta);
                  rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetTFunctor<scalar_t>;
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::R).dataView();
            auto vP = v.comp(FieldComponents::Physical::PHI).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // R component
                  rS.setProfile((-cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
                  // Phi components
                  coeff = cA*std::sin(theta);
                  rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::cos(theta);
                  rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
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
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::R).dataView();
            auto vP = v.comp(FieldComponents::Physical::THETA).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // R components
                  coeff = -cA*std::cos(theta);
                  rS.setProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::sin(theta);
                  rS.addProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
                  // Theta components
                  coeff = cA*std::sin(theta);
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::cos(theta);
                  rS.addProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPrecession::add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + preC*std::cos(alpha));
      int nR = idxFunc.dim3D();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetRFunctor>;
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::THETA).dataView();
            auto vP = v.comp(FieldComponents::Physical::PHI).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // Theta component
                  rS.addProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
                  // Phi components
                  coeff = cA*std::cos(theta);
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::sin(theta);
                  rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);

               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetTFunctor>;
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::R).dataView();
            auto vP = v.comp(FieldComponents::Physical::PHI).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // R component
                  rS.subProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
                  // Phi components
                  coeff = cA*std::sin(theta);
                  rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::cos(theta);
                  rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
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
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::R).dataView();
            auto vP = v.comp(FieldComponents::Physical::THETA).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // R components
                  coeff = cA*std::cos(theta);
                  rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::sin(theta);
                  rS.addProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
                  // Theta components
                  coeff = cA*std::sin(theta);
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::cos(theta);
                  rS.addProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);

               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPrecession::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + preC*std::cos(alpha));
      int nR = idxFunc.dim3D();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetRFunctor>;
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::THETA).dataView();
            auto vP = v.comp(FieldComponents::Physical::PHI).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // Theta component
                  rS.subProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
                  // Phi components
                  coeff = cA*std::cos(theta);
                  rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::sin(theta);
                  rS.addProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetTFunctor>;
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::R).dataView();
            auto vP = v.comp(FieldComponents::Physical::PHI).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // R component
                  rS.addProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
                  // Phi components
                  coeff = cA*std::sin(theta);
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::cos(theta);
                  rS.addProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
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
            fct_t f(cA, cB);
            Array cos_T = thGrid.array().array().cos();
            Array sin_T = thGrid.array().array().sin();
            Array cos_pt = (phGrid.array() + t).array().cos();
            Array sin_pt = (phGrid.array() + t).array().sin();
            grid_t vCosPt(const_cast<scalar_t *>(cos_pt.data()), cos_pt.size());
            grid_t vSinPt(const_cast<scalar_t *>(sin_pt.data()), sin_pt.size());
            grid_t vCosT(const_cast<scalar_t *>(cos_T.data()), cos_T.size());
            grid_t vSinT(const_cast<scalar_t *>(sin_T.data()), sin_T.size());
            Slicewise::Cpu::NoGridOp<4, fct_t, view_t, 2, 2, 0, grid_t, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            auto vT = v.comp(FieldComponents::Physical::R).dataView();
            auto vP = v.comp(FieldComponents::Physical::THETA).dataView();
            op.apply(rS.rGlobalView(), vCosT, vSinT, vCosPt, vSinPt, vT, vP, rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = idxFunc.dim2D(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(idxFunc.idx2D(iTh, iR));

                  // R components
                  coeff = cA*std::cos(theta);
                  rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::sin(theta);
                  rS.subProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
                  // Theta components
                  coeff = cA*std::sin(theta);
                  rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
                  coeff = cB*std::cos(theta);
                  rS.subProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
               }
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALPRECESSION_HPP
