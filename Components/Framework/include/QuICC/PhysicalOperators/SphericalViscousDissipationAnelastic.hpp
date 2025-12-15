/**
 * @file SphericalViscousDissipationAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

 #ifndef QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP
 #define QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP

 // System includes
 //

 // Project includes
 //
 #include "Types/Typedefs.hpp"
 #include "QuICC/Enums/FieldIds.hpp"
 #include "QuICC/VectorFields/VectorField.hpp"
 #include "QuICC/Resolutions/Resolution.hpp"
 #include "QuICC/ScalarFields/ScalarField.hpp"
 #include "QuICC/Equations/EquationParameters.hpp"
#include "DenseSM/IGenericProfile.hpp"
#include <iostream>

 namespace QuICC {

 namespace Physical {

    /**
     * @brief Implementation of the spherical coriolis term
     */
    class SphericalViscousDissipationAnelastic
    {
       public:
          /**
           * @brief Set S to viscous dissipation term
           */
          template <typename TFIELD>
          static void set(TFIELD &rS,
                          const Resolution& res,
                          const Array& r,
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);

          /**
           * @brief Add viscous dissipation term to S
           */
          template <typename TFIELD>
          static void add(TFIELD &rS,
                          const Resolution& res,
                          const Array& r,
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);

          /**
           * @brief Substract viscous dissipation term from S
           */
          template <typename TFIELD>
          static void sub(TFIELD &rS,
                          const Resolution& res,
                          const Array& r,
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);

          template <typename TFIELD>
         static void test(TFIELD &rS,
                          const Resolution& res,
                          const Array& r,
                          const Array& thGrid,    // Add theta grid
                          const Array& phGrid,    // Add phi grid
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF,
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);

          /**
           * @brief Set S to viscous dissipation term
           */
          template <typename TFIELD, typename TIDXFUNC>
          static void set(TFIELD &rS,
                          const TIDXFUNC& idxFunc,
                          const Array& nu,
                          const Array& T,
                          const Array& rho,
                          const Array& dLogRho,
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          const MHDFloat c = 1.0);

          /**
           * @brief Add viscous dissipation term to S
           */
          template <typename TFIELD, typename TIDXFUNC>
          static void add(TFIELD &rS,
                          const TIDXFUNC& idxFunc,
                          const Array& nu,
                          const Array& T,
                          const Array& rho,
                          const Array& dLogRho,
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          const MHDFloat c = 1.0);

          /**
           * @brief Substract viscous dissipation term from S
           */
          template <typename TFIELD, typename TIDXFUNC>
          static void sub(TFIELD &rS,
                          const TIDXFUNC& idxFunc,
                          const Array& nu,
                          const Array& T,
                          const Array& rho,
                          const Array& dLogRho,
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          const MHDFloat c = 1.0);

         template <typename TFIELD, typename TIDXFUNC>
         static void test(TFIELD &rS,
                          const TIDXFUNC& idxFunc,
                          const Array& r,
                          const Array& nu,
                          const Array& T,
                          const Array& rho,
                          const Array& dLogRho,
                          const Array& thGrid,    // Add theta grid
                          const Array& phGrid,    // Add phi grid
                          const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                          const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                          const MHDFloat c = 1.0);

       protected:

       private:
          /**
           * @brief Empty constructor
           */
          SphericalViscousDissipationAnelastic() = default;

          /**
           * @brief Empty destructor
           */
          ~SphericalViscousDissipationAnelastic() = default;

          /**
           * @brief helper function to calculate viscous dissipation slice
           */
         template <typename TFIELD>
          static Eigen::Matrix<MHDFloat,
                               Eigen::Dynamic,
                               Eigen::Dynamic> computeViscousSlice(const int iR,
                                                                     const int iR_,
                                                                     const MHDFloat c,
                                                                     const Datatypes::VectorField<TFIELD,
                                                                           FieldComponents::Physical::Id>& v,
                                                                     const Datatypes::TensorField<TFIELD,
                                                                           FieldComponents::Physical::Id>& Dv,
                                                                     const MHDFloat nu,
                                                                     const MHDFloat T,
                                                                     const MHDFloat rho,
                                                                     const MHDFloat dLogRho);

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

    template<typename TFIELD>
   void SphericalViscousDissipationAnelastic::set(TFIELD &rS,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto nu        = pV->evaluateLP(r, 0, 0);
      auto T         = pT->evaluateLP(r, 0, 0);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      set(rS, f, nu, T, rho, dLogRho, v, Dv, c);
   }

    template<typename TFIELD>
   void SphericalViscousDissipationAnelastic::add(TFIELD &rS,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto nu        = pV->evaluateLP(r, 0, 0);
      auto T         = pT->evaluateLP(r, 0, 0);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      add(rS, f, nu, T, rho, dLogRho, v, Dv, c);
   }

   template<typename TFIELD>
   void SphericalViscousDissipationAnelastic::sub(TFIELD &rS,
                                             const Resolution& res,
                                             const Array& r,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto nu        = pV->evaluateLP(r, 0, 0);
      auto T         = pT->evaluateLP(r, 0, 0);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      sub(rS, f, nu, T, rho, dLogRho, v, Dv, c);
   }

   template<typename TFIELD>
   void SphericalViscousDissipationAnelastic::test(TFIELD &rS,
                                             const Resolution& res,
                                             const Array& r,
                                             const Array& thGrid,    // Add theta grid
                                             const Array& phGrid,    // Add phi grid
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      IdxResFunctor f(res);
      auto nu        = pV->evaluateLP(r, 0, 0);
      auto T         = pT->evaluateLP(r, 0, 0);
      auto rho       = pF->evaluateLP(r, 0, 0);
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);
      test(rS, f, r, nu, T, rho, dLogRho, thGrid, phGrid, v, Dv, c);
   }

   template<typename TFIELD, typename TIDXFUNC>
   void SphericalViscousDissipationAnelastic::set(TFIELD &rS,
                                             const TIDXFUNC& idxFunc,
                                             const Array& nu,
                                             const Array& T,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = idxFunc.idx3D(iR);

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), rho(iR_), dLogRho(iR_));

         rS.setSlice(slice, iR);


         // Test the diagonal gradient components: OK (rms is 10^-40 or so)
         /*
         rS.setSlice(((Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::THETA).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::PHI).slice(iR).array()
                        )).matrix(), iR);
         */

         // test the curl-r: OK
         // Poloidal part is ok (r curl =0)
         /*
         rS.setSlice((v.comp(FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::PHI).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::THETA).slice(iR).array()
                        ).matrix(), iR);
         */

         // test the curl-theta: OK
         /*
         rS.setSlice(((v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::PHI).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::R).slice(iR).array()
                        )).matrix(), iR);
         */

         // test the curl-phi: OK?
         /*
         rS.addSlice(((v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::THETA).slice(iR).array()
                        )).matrix(), iR);
         */

      }
   }

   template<typename TFIELD, typename TIDXFUNC>
   void SphericalViscousDissipationAnelastic::add(TFIELD &rS,
                                             const TIDXFUNC& idxFunc,
                                             const Array& nu,
                                             const Array& T,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = idxFunc.idx3D(iR);

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), rho(iR_), dLogRho(iR_));

         rS.addSlice(slice, iR);
      }


   // *** to print the result ** //
   // std::cerr << "NL(R) = "<<rS.data()<<" \n";


   }

   template<typename TFIELD, typename TIDXFUNC>
   void SphericalViscousDissipationAnelastic::sub(TFIELD &rS,
                                             const TIDXFUNC& idxFunc,
                                             const Array& nu,
                                             const Array& T,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = idxFunc.idx3D(iR);

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), rho(iR_), dLogRho(iR_));

         rS.subSlice(slice, iR);
      }
   }

   template<typename TFIELD, typename TIDXFUNC>
   void SphericalViscousDissipationAnelastic::test(TFIELD &rS,
                                             const TIDXFUNC& idxFunc,
                                             const Array& r,
                                             const Array& nu,
                                             const Array& T,
                                             const Array& rho,
                                             const Array& dLogRho,
                                             const Array& thGrid,    // Add theta grid
                                             const Array& phGrid,    // Add phi grid
                                             const Datatypes::VectorField<TFIELD, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<TFIELD, FieldComponents::Physical::Id> &Dv,
                                             const MHDFloat c)
   {
      int nR = idxFunc.dim3D();
      int iR_;

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = idxFunc.idx3D(iR);
         int nTh = idxFunc.dim2D(iR);

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), rho(iR_), dLogRho(iR_));

         rS.addSlice(slice, iR);

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
         std::cerr << "v_r = "<<v.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";

         std::cerr << "v_rr = "<<Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_rtheta = "<<Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_rphi = "<<Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_thetar = "<<Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_thetatheta = "<<Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_thetaphi = "<<Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phir = "<<Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phitheta = "<<Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phiphi = "<<Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "Q_nu = "<< (slice/c)*T(iR_)<<" \n";
         std::cerr <<" \n";
         std::cerr << "Di*Q_nu/T = "<< slice<<" \n";
         std::cerr <<" \n";
         std::cerr << "c = "<< c<<" \n";


      }


   // *** to print the result ** //
   // std::cerr << "NL(R) = "<<rS.data()<<" \n";


   }

   // Helper function to implement the calculation of Di* Q_nu/T
   // Di is the dissipation number (passed via c);
   // Q_nu = 2 nu rho (E:E -(div(v))^2/3);
   // and E is the strain rate associate to the velocity field v = u/rho
   template <typename TFIELD>
   Eigen::Matrix<MHDFloat,
                 Eigen::Dynamic,
                 Eigen::Dynamic>SphericalViscousDissipationAnelastic::computeViscousSlice(const int iR,
                                                                                          const int iR_,
                                                                                          const MHDFloat c,
                                                                                          const Datatypes::VectorField<TFIELD,
                                                                                             FieldComponents::Physical::Id>& v,
                                                                                          const Datatypes::TensorField<TFIELD,
                                                                                             FieldComponents::Physical::Id>& Dv,
                                                                                          const MHDFloat nu,
                                                                                          const MHDFloat T,
                                                                                          const MHDFloat Rho,
                                                                                          const MHDFloat dLogRho)
{
    return (c * 2 * Rho * nu * (
                                 // E_rr^2
                                 ( -v.comp(FieldComponents::Physical::R).slice(iR).array()*dLogRho/Rho
                                   + Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::R).slice(iR).array()/Rho ).pow(2)
                                 // + E_tt^2
                                 + ( Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::THETA).slice(iR).array()/Rho ).pow(2)
                                 // + E_pp^2
                                 + ( Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::PHI).slice(iR).array()/Rho ).pow(2)
                                 // + 2* (E_rt)^2
                                 + 2*( -0.5*v.comp(FieldComponents::Physical::THETA).slice(iR).array()*dLogRho/Rho
                                      + 0.5*(  Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::THETA).slice(iR).array()
                                                +  Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::R).slice(iR).array() )/Rho ).pow(2)
                                 // + 2* (E_rp)^2
                                 + 2*( -0.5*v.comp(FieldComponents::Physical::PHI).slice(iR).array()*dLogRho/Rho
                                      + 0.5*(  Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::PHI).slice(iR).array()
                                                +  Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::R).slice(iR).array() )/Rho ).pow(2)
                                 // + 2* (E_tp)^2
                                 + 2*( 0.5*(  Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::PHI).slice(iR).array()
                                                +  Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::THETA).slice(iR).array() )/Rho ).pow(2)
                                 // -(1/3)div(v)
                                 - (1.0/3.0) * ( -v.comp(FieldComponents::Physical::R).slice(iR).array()*dLogRho/Rho ).pow(2)
                              ) / T).matrix();
}
 } // namespace Physical
 } // namespace QuICC

 #endif // QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP

