/**
 * @file SphericalOhmicDissipationAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALOHMICDISSIPATIONANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALOHMICDISSIPATIONANELASTIC_HPP


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
   class SphericalOhmicDissipationAnelastic
   {
      public:
         /**
          * @brief Add Ohmic dissipation term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS,
                         const Resolution& res,
                         const Array& r,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief test Ohmic dissipation term to S
          */
         template <typename TFIELD>
         static void test(TFIELD &rS,
                         const Resolution& res,
                         const Array& r,
                         const Array& thGrid,    // Add theta grid
                         const Array& phGrid,    // Add phi grid
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Ohmic dissipation term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS,
                         const TIDXFUNC& idxFunc,
                         const Array& eta,
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

         /**
          * @brief test Ohmic dissipation term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void test(TFIELD &rS,
                         const TIDXFUNC& idxFunc,
                         const Array& r,
                         const Array& eta,
                         const Array& thGrid,    // Add theta grid
                         const Array& phGrid,    // Add phi grid
                         const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalOhmicDissipationAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalOhmicDissipationAnelastic() = default;

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
            QUICC_CUDA_HOSTDEV T operator()(T gr, T ui, T uj, T uk)
            {
               return _scaling * gr * (ui * ui  + uj * uj + uk * uk);
            }
         };

   };

   template <typename TFIELD>
   void SphericalOhmicDissipationAnelastic::add(TFIELD &rS,
                                                const Resolution& res,
                                                const Array& r,
                                                const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                                std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for magnetic diffusivity, Eta
                                                const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                                const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto eta = pF->evaluateLP(r, 0, 0);
      add(rS, f, eta, w, c);
   }

   template <typename TFIELD>
   void SphericalOhmicDissipationAnelastic::test(TFIELD &rS,
                                                const Resolution& res,
                                                const Array& r,
                                                const Array& thGrid,    // Add theta grid
                                                const Array& phGrid,    // Add phi grid
                                                const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                                std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for magnetic diffusivity, Eta
                                                const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                                const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto eta = pF->evaluateLP(r, 0, 0);
      test(rS, f, r, eta, thGrid, phGrid, w, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalOhmicDissipationAnelastic::add(TFIELD &rS,
                                                const TIDXFUNC& idxFunc,
                                                const Array& eta,
                                                const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                                const MHDFloat c)
   {
      using scalar_t = typename TFIELD::PointType;
      if constexpr(std::is_same_v<TFIELD, Datatypes::ViewScalarField<scalar_t>>)
      {
         using view_t = typename Datatypes::ViewScalarField<scalar_t>::ViewStorageType;
         using grid_t = View::ViewBase<double>;
         using fct_t = details::AddTmplFunctor<scalar_t, SetFunctor>;
         fct_t f(c);
         grid_t vEta(const_cast<scalar_t *>(eta.data()), eta.size());
         Slicewise::Cpu::NoGridOp<2, fct_t, view_t, 1, 0, 0, grid_t, view_t, view_t, view_t, view_t> op(f);
         auto vR = w.comp(FieldComponents::Physical::R).dataView();
         auto vT = w.comp(FieldComponents::Physical::THETA).dataView();
         auto vP = w.comp(FieldComponents::Physical::PHI).dataView();
         op.apply(rS.rGlobalView(), vEta, vR, vT, vP, rS.dataView());
      }
      else
      {
         int nR = idxFunc.dim3D();
         int iR_;

         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = idxFunc.idx3D(iR);

            // Boussinesq part (not vanishing for dLogEta =0)
            rS.addSlice(c*eta(iR_)*(   w.comp(FieldComponents::Physical::R).slice(iR).array()
                     * w.comp(FieldComponents::Physical::R).slice(iR).array()
                     ).matrix(), iR);

            rS.addSlice(c*eta(iR_)*(   w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                     * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                     ).matrix(), iR);

            rS.addSlice(c*eta(iR_)*(   w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                     * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                     ).matrix(), iR);

         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalOhmicDissipationAnelastic::test(TFIELD &rS,
                                                const TIDXFUNC& idxFunc,
                                                const Array& r,
                                                const Array& eta,
                                                const Array& thGrid,    // Add theta grid
                                                const Array& phGrid,    // Add phi grid
                                                const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &w,
                                                const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = idxFunc.idx3D(iR);
            int nTh = idxFunc.dim2D(iR);

            // Boussinesq part (not vanishing for dLogEta =0)
            rS.setSlice(c*eta(iR_)*(   w.comp(FieldComponents::Physical::R).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    ).matrix(), iR);

            rS.addSlice(c*eta(iR_)*(   w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    ).matrix(), iR);

            rS.addSlice(c*eta(iR_)*(   w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    ).matrix(), iR);

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
            std::cerr << "eta(iR) =  ("<<eta(iR)<<")"<<" \n";
            std::cerr <<" \n";

            std::cerr << "j_r = "<<w.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "j_theta = "<<w.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "j_phi = "<<w.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "Q_j = "<< eta(iR_)*(   w.comp(FieldComponents::Physical::R).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array() +
                                    w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array() +
                              w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()).matrix() <<" \n";
            std::cerr <<" \n";
            std::cerr << "c = "<< c<<" \n";
         }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALOHMICDISSIPATIONANELASTIC_HPP
