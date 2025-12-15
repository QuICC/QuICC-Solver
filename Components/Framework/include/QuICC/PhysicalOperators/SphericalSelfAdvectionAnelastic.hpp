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
   };

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

      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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

      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP
