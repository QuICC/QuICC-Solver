/** 
 * @file LegendreBase.cpp
 * @brief Source of the implementation of the Legendre polynomial
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Legendre/LegendreBase.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Legendre {

   LegendreBase::LegendreBase()
   {
   }

   LegendreBase::~LegendreBase()
   {
   }

   void LegendreBase::Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerN norm)
   {
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn);

      ipn.array() = cs(1)*(cs(0)*ipn_2.array() + (igrid.array()*ipn_1.array()).array());
   }
   
   void LegendreBase::P0(Eigen::Ref<internal::Matrix> ip0, Normalizer norm)
   {
      internal::Array cs = norm();

      ip0.setConstant(cs(0));
   }

   void LegendreBase::P1(Eigen::Ref<internal::Matrix> ip1, const internal::Array& igrid, Normalizer norm)
   {
      internal::Array cs = norm();

      ip1.array() = (cs(0)*igrid);
   }

   void LegendreBase::dP1(Eigen::Ref<internal::Matrix> idp1, const Eigen::Ref<const internal::Matrix>& ipn_1, Normalizer norm)
   {
      internal::Array cs = norm();
 
      idp1.array() = cs(0)*ipn_1.array();
   }

   void LegendreBase::dPn(Eigen::Ref<internal::Matrix> idpn, const int n, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& idpn_1, const internal::Array& igrid, NormalizerN norm)
   {
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn);
 
      idpn.array() = cs(0)*ipn_1.array() + cs(1)*(igrid.array()*idpn_1.array()).array();
   }

   //
   // General polynomial normalizer
   //
   LegendreBase::Normalizer LegendreBase::normP0()
   {
      #if defined QUICC_LEGENDRE_NORM_NATURAL
         return &LegendreBase::naturalP0;
      #elif defined QUICC_LEGENDRE_NORM_UNITY
         return &LegendreBase::unitP0;
      #endif //defined QUICC_LEGENDRE_NORM_NATURAL
   }

   LegendreBase::Normalizer LegendreBase::normP1()
   {
      #if defined QUICC_LEGENDRE_NORM_NATURAL
         return &LegendreBase::naturalP1;
      #elif defined QUICC_LEGENDRE_NORM_UNITY
         return &LegendreBase::unitP1;
      #endif //defined QUICC_LEGENDRE_NORM_NATURAL
   }

   LegendreBase::NormalizerN LegendreBase::normPn()
   {
      #if defined QUICC_LEGENDRE_NORM_NATURAL
         return &LegendreBase::naturalPn;
      #elif defined QUICC_LEGENDRE_NORM_UNITY
         return &LegendreBase::unitPn;
      #endif //defined QUICC_LEGENDRE_NORM_NATURAL
   }

   LegendreBase::Normalizer LegendreBase::normdP1()
   {
      #if defined QUICC_LEGENDRE_NORM_NATURAL
         return &LegendreBase::naturaldP1;
      #elif defined QUICC_LEGENDRE_NORM_UNITY
         return &LegendreBase::unitdP1;
      #endif //defined QUICC_LEGENDRE_NORM_NATURAL
   }

   LegendreBase::NormalizerN LegendreBase::normdPn()
   {
      #if defined QUICC_LEGENDRE_NORM_NATURAL
         return &LegendreBase::naturaldPn;
      #elif defined QUICC_LEGENDRE_NORM_UNITY
         return &LegendreBase::unitdPn;
      #endif //defined QUICC_LEGENDRE_NORM_NATURAL
   }

   internal::Array LegendreBase::unitP0()
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(0.5));

      return cs;
   }

   internal::Array LegendreBase::unitP1()
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(3.0)/MHD_MP(2.0));

      return cs;
   }

   internal::Array LegendreBase::unitPn(const internal::MHDFloat dn)
   {
      internal::Array cs(2);

      cs(0) = -precision::sqrt((MHD_MP(2.0)*dn - MHD_MP(1.0))/(MHD_MP(2.0)*dn - MHD_MP(3.0)))*(dn - MHD_MP(1.0))/(MHD_MP(2.0)*dn - MHD_MP(1.0));
      cs(1) = precision::sqrt((MHD_MP(2.0)*dn + MHD_MP(1.0))/(MHD_MP(2.0)*dn - MHD_MP(1.0)))*(MHD_MP(2.0)*dn-MHD_MP(1.0))/dn;

      return cs;
   }

   internal::Array LegendreBase::unitdP1()
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(3.0));

      return cs;
   }

   internal::Array LegendreBase::unitdPn(const internal::MHDFloat dn)
   {
      internal::Array cs(2);

      cs(0) = dn*precision::sqrt((MHD_MP(2.0)*dn + MHD_MP(1.0))/(MHD_MP(2.0)*dn - MHD_MP(1.0)));

      cs(1) = precision::sqrt((MHD_MP(2.0)*dn + MHD_MP(1.0))/(MHD_MP(2.0)*dn - MHD_MP(1.0)));

      return cs;
   }

   internal::Array LegendreBase::naturalP0()
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   internal::Array LegendreBase::naturalP1()
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   internal::Array LegendreBase::naturalPn(const internal::MHDFloat dn)
   {
      internal::Array cs(2);

      cs(0) = -(dn - MHD_MP(1.0))/(MHD_MP(2.0)*dn - MHD_MP(1.0));
      cs(1) = (MHD_MP(2.0)*dn - MHD_MP(1.0))/dn;

      return cs;
   }

   internal::Array LegendreBase::naturaldP1()
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   internal::Array LegendreBase::naturaldPn(const internal::MHDFloat dn)
   {
      internal::Array cs(2);

      cs(0) = dn;

      cs(1) = MHD_MP(1.0);

      return cs;
   }

}
}
}
