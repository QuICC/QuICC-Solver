/**
 * @file SphericalSelfAdvectionAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "DenseSM/IGenericProfile.hpp"
#include "ViewOps/Slicewise/NoGridOp.hpp"
#include "QuICC/PhysicalOperators/details/FunctorHelpers.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalSelfAdvectionAnelastic
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
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD>
         static void sub(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Set S to Coriolis term
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void set(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Array& dLogRho,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Array& dLogRho,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void sub(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Array& dLogRho,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalSelfAdvectionAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalSelfAdvectionAnelastic() = default;

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
         template <class T = double> struct SetRFunctor
         {
            /// @brief non dimensional scaling for transport term
            T _scaling;

            /// @brief ctor
            /// @param scaling
            SetRFunctor(T scaling) : _scaling(scaling){};

            /// @brief deleted default constructor
            SetRFunctor() = delete;

            /// @brief dtor
            ~SetRFunctor() = default;

            /// @brief Dot product
            /// @param g
            /// @param ui
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T rho, T dLog, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return _scaling * (uk * vj - uj * vk + dLog * (uj * uj + uk * uk)) / (rho * rho);
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
            /// @param gc
            /// @param gs
            /// @param ui
            /// @param uj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T rho, T dLog, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return _scaling * (ui * vk - uk * vi - dLog * (ui * uj)) / (rho * rho);
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
            /// @param gc
            /// @param gs
            /// @param ui
            /// @param uj
            /// @return
            QUICC_CUDA_HOSTDEV T operator()(T rho, T dLog, T ui, T uj, T uk, T vi, T vj, T vk)
            {
               return _scaling * (uj * vi - ui * vj - dLog * (ui * uk )) / (rho * rho);
            }
         };

         template <typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <typename TFIELD> void SphericalSelfAdvectionAnelastic::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f)
   {
      for(auto&& [k, flag]: f.enabled())
      {
         vs.push_back(f.comp(k).dataView());
      }
   }

   template <typename TFIELD>
   void SphericalSelfAdvectionAnelastic::set(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // intended for derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      set(rS, compId, f, rho, dLogRho, v, w, c);
   }

   template <typename TFIELD>
   void SphericalSelfAdvectionAnelastic::add(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // intended for derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      add(rS, compId, f, rho, dLogRho, v, w, c);
   }

   template <typename TFIELD>
   void SphericalSelfAdvectionAnelastic::sub(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // intended for derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      sub(rS, compId, f, rho, dLogRho, v, w, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalSelfAdvectionAnelastic::set(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const TIDXFUNC& idxFunc,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetRFunctor<scalar_t>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5]);
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.setSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

                  rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              )/ rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.setSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

                  rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              )/ rho(iR_) / rho(iR_)).matrix(),  iR);

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
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5]);
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.setSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.setSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
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
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5]);
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.setSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.setSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalSelfAdvectionAnelastic::add(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const TIDXFUNC& idxFunc,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            //using fct_t = AddRFunctor<scalar_t>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetRFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

                  rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              )/ rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.addSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

                  rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              )/ rho(iR_) / rho(iR_)).matrix(),  iR);

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
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.addSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
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
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.addSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.subSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);
               }
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalSelfAdvectionAnelastic::sub(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const TIDXFUNC& idxFunc,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::SubTmplFunctor<scalar_t, SetRFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

                  rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              )/ rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.subSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.addSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

                  rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              )/ rho(iR_) / rho(iR_)).matrix(),  iR);
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
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.subSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.addSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

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
            grid_t vDLog(const_cast<scalar_t *>(dLogRho.data()), dLogRho.size());
            std::vector<view_t> vs;
            collectViews(vs, v);
            collectViews(vs, w);
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 2, grid_t, grid_t, view_t, view_t, view_t, view_t, view_t, view_t, view_t> op(f);
            op.apply(rS.rGlobalView(), vRho, vDLog, vs[0], vs[1], vs[2], vs[3], vs[4], vs[5], rS.dataView());
         }
         else
         {
            if(c != 1.0)
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            } else
            {
               for(int iR = 0; iR < nR; ++iR)
               {
                  iR_ = idxFunc.idx3D(iR);

                  // Boussinesq part (not vanishing for dLogRho =0)
                  rS.subSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           * w.comp(FieldComponents::Physical::R).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  rS.addSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array()
                           * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                           / rho(iR_) / rho(iR_)).matrix(), iR);

                  // Anelastic part (vanishing for dLogRho =0)
                  rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array()
                              * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                              ) / rho(iR_) / rho(iR_)).matrix(),  iR);

               }
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP
