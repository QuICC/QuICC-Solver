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
   };

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


      } else if(compId == FieldComponents::Physical::THETA)
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

      } else if(compId == FieldComponents::Physical::PHI)
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


      } else if(compId == FieldComponents::Physical::THETA)
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

      } else if(compId == FieldComponents::Physical::PHI)
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
