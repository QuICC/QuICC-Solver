/*
 * @file SphericalPoincare.hpp
 * @brief Implementation of the spherical poincare term
 */

#ifndef QUICC_PHYSICAL_SPHERICALPOINCARE_HPP
#define QUICC_PHYSICAL_SPHERICALPOINCARE_HPP

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
    * @brief Implementation of the spherical poincare term
    */
   class SphericalPoincare
   {
      public:
         /**
          * @brief Set S to Poincare term
          */
         template <typename TFIELD>
         static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Add Poincare term to S
          */
         template <typename TFIELD>
         static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Substract Poincare term from S
          */
         template <typename TFIELD>
         static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Set S to Poincare term
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Add Poincare term to S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Substract Poincare term from S
          */
         template <typename TFIELD, typename TIDXFUNC>
         static void sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPoincare() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalPoincare() = default;

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
   void SphericalPoincare::set(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      IdxResFunctor f(res);
      set(rS, compId, f, rGrid, thGrid, phGrid, c);
   }

   template <typename TFIELD>
   void SphericalPoincare::add(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      IdxResFunctor f(res);
      add(rS, compId, f, rGrid, thGrid, phGrid, c);
   }

   template <typename TFIELD>
   void SphericalPoincare::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      IdxResFunctor f(res);
      sub(rS, compId, f, rGrid, thGrid, phGrid, c);
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPoincare::set(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         rS.setZeros();
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = -csa*rGrid(idxFunc.idx3D(iR));
            nTh = idxFunc.dim2D(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               rS.setProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(idxFunc.idx3D(iR));
            nTh = idxFunc.dim2D(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = idxFunc.idx2D(iTh, iR);

               rS.setProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPoincare::add(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         //
         // Zero
         //
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(idxFunc.idx3D(iR));
            nTh = idxFunc.dim2D(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(idxFunc.idx3D(iR));
            nTh = idxFunc.dim2D(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = idxFunc.idx2D(iTh, iR);

               rS.addProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
            }
         }
      }
   }

   template <typename TFIELD, typename TIDXFUNC>
   void SphericalPoincare::sub(TFIELD &rS, FieldComponents::Physical::Id compId, const TIDXFUNC& idxFunc, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = idxFunc.dim3D();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         //
         // Zero
         //
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(idxFunc.idx3D(iR));
            nTh = idxFunc.dim2D(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(idxFunc.idx3D(iR));
            nTh = idxFunc.dim2D(iR);
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = idxFunc.idx2D(iTh, iR);

               rS.subProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
            }
         }
      }
   }
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_SPHERICALPOINCARE_HPP
