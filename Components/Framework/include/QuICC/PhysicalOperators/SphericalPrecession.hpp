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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
      } else if(compId == FieldComponents::Physical::THETA)
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
      } else if(compId == FieldComponents::Physical::PHI)
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
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALPRECESSION_HPP
