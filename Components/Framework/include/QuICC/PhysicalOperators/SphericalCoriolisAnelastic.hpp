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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALCORIOLISANELASTIC_HPP
