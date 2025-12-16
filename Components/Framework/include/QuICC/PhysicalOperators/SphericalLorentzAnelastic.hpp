/**
 * @file SphericalLorentzAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALLORENTZANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALLORENTZANELASTIC_HPP

// System includes
//
#include <iostream>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "DenseSM/IGenericProfile.hpp"
#include "QuICC/PhysicalOperators/details/FunctorHelpers.hpp"
#include "ViewOps/Slicewise/Cpu/NoGridOp.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalLorentzAnelastic
   {
      public:
         /**
          * @brief Set Lorentz term to S
          */
         template <typename TFIELD>
         static void set(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Lorentz term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Lorentz term to S
          */
         template <typename TFIELD>
         static void test(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const Resolution& res,
                         const Array& r,
                         const Array& thGrid,    // Add theta grid
                         const Array& phGrid,    // Add phi grid
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Set Lorentz term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void set(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Lorentz term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& rho,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Lorentz term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void test(TFIELD &rS,
                         FieldComponents::Physical::Id compId,
                         const TIDXFUNC& idxFunc,
                         const Array& r,
                         const Array& rho,
                         const Array& thGrid,    // Add theta grid
                         const Array& phGrid,    // Add phi grid
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalLorentzAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalLorentzAnelastic() = default;

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
            QUICC_CUDA_HOSTDEV T operator()(T gr, T ui, T uj, T vi, T vj)
            {
               return _scaling * (ui * vj  - uj * vi) / gr;
            }
         };

         template <FieldComponents::Physical::Id C1, FieldComponents::Physical::Id C2, typename TFIELD>
         static void collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD& f);
   };

   template <FieldComponents::Physical::Id C1, FieldComponents::Physical::Id C2, typename TFIELD>
      void SphericalLorentzAnelastic::collectViews(std::vector<typename TFIELD::ScalarFieldType::ViewStorageType>& vs, const TFIELD &f)
      {
         vs.emplace_back(f.comp(C1).dataView());
         vs.emplace_back(f.comp(C2).dataView());
      }

   template <typename TFIELD>
   void SphericalLorentzAnelastic::set(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, rho
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho = pF->evaluateLP(r, 0, 0);
      set(rS, compId, f, rho, v, w, c);
   }

   template <typename TFIELD>
   void SphericalLorentzAnelastic::add(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, rho
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho = pF->evaluateLP(r, 0, 0);
      add(rS, compId, f, rho, v, w, c);
   }

   template <typename TFIELD>
   void SphericalLorentzAnelastic::test(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const Resolution& res,
                                             const Array& r,
                                             const Array& thGrid,    // Add theta grid
                                             const Array& phGrid,    // Add phi grid
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, rho
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto rho = pF->evaluateLP(r, 0, 0);
      add(rS, compId, f, r, rho, thGrid, phGrid, v, w, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalLorentzAnelastic::set(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const TIDXFUNC& idxFunc,
                                             const Array& rho,
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
            using fct_t = SetFunctor<scalar_t>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t> op(f);
            std::vector<view_t> vs;
            vs.reserve(4);
            collectViews<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>(vs, v);
            collectViews<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>(vs, w);
            op.apply(rS.rGlobalView(), vRho, vs[0], vs[1], vs[2], vs[3]);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);

               // Boussinesq part (not vanishing for dLogrho =0)
               rS.setSlice(-c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetFunctor<scalar_t>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t> op(f);
            std::vector<view_t> vs;
            vs.reserve(4);
            collectViews<FieldComponents::Physical::PHI,FieldComponents::Physical::R>(vs, v);
            collectViews<FieldComponents::Physical::PHI,FieldComponents::Physical::R>(vs, w);
            op.apply(rS.rGlobalView(), vRho, vs[0], vs[1], vs[2], vs[3]);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);

               // Boussinesq part (not vanishing for dLogrho =0)
               rS.setSlice(-c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                        * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        * w.comp(FieldComponents::Physical::R).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

            }
         }

      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = SetFunctor<scalar_t>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t> op(f);
            std::vector<view_t> vs;
            vs.reserve(4);
            collectViews<FieldComponents::Physical::R,FieldComponents::Physical::THETA>(vs, v);
            collectViews<FieldComponents::Physical::R,FieldComponents::Physical::THETA>(vs, w);
            op.apply(rS.rGlobalView(), vRho, vs[0], vs[1], vs[2], vs[3]);
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);

               // Boussinesq part (not vanishing for dLogrho =0)
               rS.setSlice(-c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        * w.comp(FieldComponents::Physical::R).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                        * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        / rho(iR_)).matrix(), iR);
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalLorentzAnelastic::add(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const TIDXFUNC& idxFunc,
                                             const Array& rho,
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
            using fct_t = details::AddTmplFunctor<scalar_t, SetFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t, view_t> op(f);
            std::vector<view_t> vs;
            vs.reserve(4);
            collectViews<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>(vs, v);
            collectViews<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>(vs, w);
            op.apply(rS.rGlobalView(), vRho, vs[0], vs[1], vs[2], vs[3], rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);

               // Boussinesq part (not vanishing for dLogrho =0)
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t, view_t> op(f);
            std::vector<view_t> vs;
            vs.reserve(4);
            collectViews<FieldComponents::Physical::PHI,FieldComponents::Physical::R>(vs, v);
            collectViews<FieldComponents::Physical::PHI,FieldComponents::Physical::R>(vs, w);
            op.apply(rS.rGlobalView(), vRho, vs[0], vs[1], vs[2], vs[3], rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);

               // Boussinesq part (not vanishing for dLogrho =0)
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                        * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        * w.comp(FieldComponents::Physical::R).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

            }
         }

      } else if(compId == FieldComponents::Physical::PHI)
      {
         using scalar_t = typename TFIELD::PointType;
         if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
         {
            using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
            using grid_t = View::ViewBase<double>;
            using fct_t = details::AddTmplFunctor<scalar_t, SetFunctor>;
            fct_t f(c);
            grid_t vRho(const_cast<scalar_t *>(rho.data()), rho.size());
            Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t, view_t> op(f);
            std::vector<view_t> vs;
            vs.reserve(4);
            collectViews<FieldComponents::Physical::R,FieldComponents::Physical::THETA>(vs, v);
            collectViews<FieldComponents::Physical::R,FieldComponents::Physical::THETA>(vs, w);
            op.apply(rS.rGlobalView(), vRho, vs[0], vs[1], vs[2], vs[3], rS.dataView());
         }
         else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = idxFunc.idx3D(iR);

               // Boussinesq part (not vanishing for dLogrho =0)
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        * w.comp(FieldComponents::Physical::R).slice(iR).array()
                        / rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                        * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        / rho(iR_)).matrix(), iR);
            }

         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalLorentzAnelastic::test(TFIELD &rS,
                                             FieldComponents::Physical::Id compId,
                                             const TIDXFUNC& idxFunc,
                                             const Array& r,
                                             const Array& rho,
                                             const Array& thGrid,    // Add theta grid
                                             const Array& phGrid,    // Add phi grid
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = idxFunc.idx3D(iR);
            int nTh = idxFunc.dim2D(iR);

            // Boussinesq part (not vanishing for dLogrho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / rho(iR_)).matrix(), iR);

            // to test the implementation
            std::cerr << "(iR, iR_) = ("<<iR<<","<<iR_<<")"<<" \n";
            std::cerr << "(r(iR), r(iR_)) = ("<<r(iR)<<","<<r(iR_)<<")"<<" \n";
            // Print theta and phi coordinates
            std::cerr << " theta = \n";

            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               int iTh_ = idxFunc.idx2D(iTh, iR);
               MHDFloat theta = thGrid(iTh_);

               std::cerr << theta << " ";

            }

            // Print phi values - use grid size directly since phi is typically uniform
            int nPh = phGrid.size();
            std::cerr << "\n phi = \n";
            for(int iPh = 0; iPh < nPh; ++iPh)
            {
               MHDFloat phi = phGrid(iPh);
               std::cerr <<phi << " ";
            }

            std::cerr << "\n";
            // for density_type=0, this is r
            std::cerr << "rho(iR) =  ("<<rho(iR)<<")"<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "w_theta = "<<w.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "w_phi = "<<w.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "rnLcomp = "<< -c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / rho(iR_)).matrix() + c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / rho(iR_)).matrix() <<" \n";


         }


      } else if(compId == FieldComponents::Physical::THETA)
      {

         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = idxFunc.idx3D(iR);
            int nTh = idxFunc.dim2D(iR);

            // Boussinesq part (not vanishing for dLogrho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / rho(iR_)).matrix(), iR);

            // to test the implementation
            std::cerr << "(iR, iR_) = ("<<iR<<","<<iR_<<")"<<" \n";
            std::cerr << "(r(iR), r(iR_)) = ("<<r(iR)<<","<<r(iR_)<<")"<<" \n";
            // Print theta and phi coordinates
            std::cerr << " theta = \n";

            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               int iTh_ = idxFunc.idx2D(iTh, iR);
               MHDFloat theta = thGrid(iTh_);

               std::cerr << theta << " ";

            }

            // Print phi values - use grid size directly since phi is typically uniform
            int nPh = phGrid.size();
            std::cerr << "\n phi = \n";
            for(int iPh = 0; iPh < nPh; ++iPh)
            {
               MHDFloat phi = phGrid(iPh);
               std::cerr <<phi << " ";
            }

            std::cerr << "\n";
            // for density_type=0, this is r
            std::cerr << "rho(iR) =  ("<<rho(iR)<<")"<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "w_theta = "<<w.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "w_phi = "<<w.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "rnLcomp = "<< -c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / rho(iR_)).matrix() + c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / rho(iR_)).matrix() <<" \n";



         }

      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = idxFunc.idx3D(iR);
            int nTh = idxFunc.dim2D(iR);

            // Boussinesq part (not vanishing for dLogrho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / rho(iR_)).matrix(), iR);

            // to test the implementation
            std::cerr << "(iR, iR_) = ("<<iR<<","<<iR_<<")"<<" \n";
            std::cerr << "(r(iR), r(iR_)) = ("<<r(iR)<<","<<r(iR_)<<")"<<" \n";
            // Print theta and phi coordinates
            std::cerr << " theta = \n";

            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               int iTh_ = idxFunc.idx2D(iTh, iR);
               MHDFloat theta = thGrid(iTh_);

               std::cerr << theta << " ";

            }

            // Print phi values - use grid size directly since phi is typically uniform
            int nPh = phGrid.size();
            std::cerr << "\n phi = \n";
            for(int iPh = 0; iPh < nPh; ++iPh)
            {
               MHDFloat phi = phGrid(iPh);
               std::cerr <<phi << " ";
            }

            std::cerr << "\n";
            // for density_type=0, this is r
            std::cerr << "rho(iR) =  ("<<rho(iR)<<")"<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "w_theta = "<<w.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "w_phi = "<<w.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "rnLcomp = "<< -c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / rho(iR_)).matrix() + c*(   v.comp(FieldComponents::Physical::R).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / rho(iR_)).matrix() <<" \n";


         }

      }

   }
} // namespace Physical
} // namespace QuIC

#endif // QUICC_PHYSICAL_SPHERICALLORENTZANELASTIC_HPP
