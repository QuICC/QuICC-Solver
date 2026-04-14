/**
 * @file SphericalCoriolisAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "DenseSM/IGenericProfile.hpp"
#include "QuICC/PhysicalOperators/details/FunctorHelpers.hpp"
#include "ViewOps/Slicewise/Cpu/NoGridOp.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalCoriolisAnelastic
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD>
         static void set(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Array& cosTheta,
                         const Array& sinTheta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Array& cosTheta,
                         const Array& sinTheta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD>
         static void sub(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Array& cosTheta,
                         const Array& sinTheta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const MHDFloat c = 1.0);

         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void set(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Array& cosTheta,
                         const Array& sinTheta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Array& cosTheta,
                         const Array& sinTheta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void sub(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Array& cosTheta,
                         const Array& sinTheta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalCoriolisAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalCoriolisAnelastic() = default;

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
         template <class T = double> struct SetRTFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetRTFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetRTFunctor() = delete;

            /// @brief dtor
            ~SetRTFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gr, T gt, T ui)
            {
               return -_scaling * (gt * ui / gr);
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

            ///@brief deleted default constructor
            SetPFunctor() = delete;

            /// @brief dtor
            ~SetPFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T gr, T gc, T gs, T ui, T uj)
            {
               return _scaling * (gs * ui + gc * uj) / gr;
            }
         };
   };

   template <typename TFIELD>
   void SphericalCoriolisAnelastic::set(TFIELD &rS,
                                        FieldComponents::Physical::Id compId,
                                        const Resolution& res,
                                        const Array& r,
                                        const Array& cosTheta,
                                        const Array& sinTheta,
                                        const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                                        const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho = pF->evaluateLP(r, 0, 0);
      set(rS, compId, f, rho, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD>
   void SphericalCoriolisAnelastic::add(TFIELD &rS,
                                        FieldComponents::Physical::Id compId,
                                        const Resolution& res,
                                        const Array& r,
                                        const Array& cosTheta,
                                        const Array& sinTheta,
                                        const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                                        const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho = pF->evaluateLP(r, 0, 0);
      add(rS, compId, f, rho, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD>
   void SphericalCoriolisAnelastic::sub(TFIELD &rS,
                                        FieldComponents::Physical::Id compId,
                                        const Resolution& res,
                                        const Array& r,
                                        const Array& cosTheta,
                                        const Array& sinTheta,
                                        const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                                        const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho = pF->evaluateLP(r, 0, 0);
      sub(rS, compId, f, rho, cosTheta, sinTheta, v, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalCoriolisAnelastic::set(TFIELD &rS,
                                        FieldComponents::Physical::Id compId,
                                        const TIDXFUNC& idxFunc,
                                        const Array& rho,
                                        const Array& cosTheta,
                                        const Array& sinTheta,
                                        const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                        const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vSin, v.comp(FieldComponents::Physical::PHI).dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)/rho(iR_), iTh, iR);
                  }
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
            using fct_t = SetRTFunctor<scalar_t>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vCos, v.comp(FieldComponents::Physical::PHI).dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(-v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)/rho(iR_), iTh, iR);
                  }
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
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 2, 0, grid_t, grid_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vCos, vSin, v.comp(FieldComponents::Physical::R).dataView(), v.comp(FieldComponents::Physical::THETA).dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);
                     rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.setProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)/rho(iR_), iTh, iR);
                     rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)/rho(iR_), iTh, iR);
                  }
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalCoriolisAnelastic::add(TFIELD &rS,
                                        FieldComponents::Physical::Id compId,
                                        const TIDXFUNC& idxFunc,
                                        const Array& rho,
                                        const Array& cosTheta,
                                        const Array& sinTheta,
                                        const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                        const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetRTFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vSin, v.comp(FieldComponents::Physical::PHI).dataView(), rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);

                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)/rho(iR_), iTh, iR);
                  }
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
            using fct_t = details::AddTmplFunctor<scalar_t, SetRTFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vCos, v.comp(FieldComponents::Physical::PHI).dataView(), rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)/rho(iR_), iTh, iR);
                  }
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
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 2, 0, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vCos, vSin, v.comp(FieldComponents::Physical::R).dataView(), v.comp(FieldComponents::Physical::THETA).dataView(), rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);
                     rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)/rho(iR_), iTh, iR);
                     rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)/rho(iR_), iTh, iR);
                  }
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalCoriolisAnelastic::sub(TFIELD &rS,
                                        FieldComponents::Physical::Id compId,
                                        const TIDXFUNC& idxFunc,
                                        const Array& rho,
                                        const Array& cosTheta,
                                        const Array& sinTheta,
                                        const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                        const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetRTFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vSin, v.comp(FieldComponents::Physical::PHI).dataView(), rS.dataView());
         }
         else
         {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               nTh = idxFunc.dim2D(iR);

               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.addProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);
               }
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);
               nTh = idxFunc.dim2D(iR);

               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = idxFunc.idx2D(iTh, iR);

                  rS.addProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)/rho(iR_), iTh, iR);
               }
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
            using fct_t = details::SubTmplFunctor<scalar_t, SetRTFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 1, 0, grid_t, grid_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vCos, v.comp(FieldComponents::Physical::PHI).dataView(), rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.addProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)/rho(iR_), iTh, iR);
                  }
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
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vCos(const_cast<scalar_t *>(cosTheta.data()), cosTheta.size());
            grid_t vSin(const_cast<scalar_t *>(sinTheta.data()), sinTheta.size());
            Slicewise::Cpu::NoGridOp<3, fct_t, view_t, 1, 2, 0, grid_t, grid_t, grid_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vCos, vSin, v.comp(FieldComponents::Physical::R).dataView(), v.comp(FieldComponents::Physical::THETA).dataView(), rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);
                     rS.subProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_))/rho(iR_), iTh, iR);
                  }
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);
                  nTh = idxFunc.dim2D(iR);

                  for(int iTh = 0; iTh < nTh; ++iTh)
                  {
                     iTh_ = idxFunc.idx2D(iTh, iR);

                     rS.subProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)/rho(iR_), iTh, iR);
                     rS.subProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)/rho(iR_), iTh, iR);
                  }
               }
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP
