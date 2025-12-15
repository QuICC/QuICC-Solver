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
