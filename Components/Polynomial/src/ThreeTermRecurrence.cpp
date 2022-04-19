/**
 * @file ThreeTermRecurrence.cpp
 * @brief Source of the general three term recurrence implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

   void ThreeTermRecurrence::Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerNC norm)
   {
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn);

      ipn.array() = cs(3)*(cs(0)*ipn_2.array() + (cs(1)*igrid.array() + cs(2))*ipn_1.array());
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<internal::Matrix> ip1, const Eigen::Ref<const internal::Matrix>& ip0, const internal::Array& igrid, NormalizerC norm)
   {
      internal::Array cs = norm();

      ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1))*ip0.array();
   }

   void ThreeTermRecurrence::P0(Eigen::Ref<internal::Matrix> ip0, const internal::Array& igrid, NormalizerC norm)
   {
      internal::Array cs = norm();

      ip0.setConstant(cs(0));
   }

   void ThreeTermRecurrence::Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const internal::MHDFloat alpha, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerNA norm)
   {
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn, alpha);

      ipn.array() = cs(3)*(cs(0)*ipn_2.array() + (cs(1)*igrid.array() + cs(2))*ipn_1.array());
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<internal::Matrix> ip1, const internal::MHDFloat alpha, const Eigen::Ref<const internal::Matrix>& ip0, const internal::Array& igrid, NormalizerA norm)
   {
      internal::Array cs = norm(alpha);

      ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1))*ip0.array();
   }

   void ThreeTermRecurrence::P0(Eigen::Ref<internal::Matrix> ip0, const internal::MHDFloat alpha, const internal::Array& igrid, NormalizerA norm)
   {
      internal::Array cs = norm(alpha);

      ip0.setConstant(cs(0));
   }

   void ThreeTermRecurrence::Pn(Eigen::Ref<internal::Matrix> ipn, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipn_1, const Eigen::Ref<const internal::Matrix>& ipn_2, const internal::Array& igrid, NormalizerNAB norm)
   {
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn, alpha, beta);

      if(cs(2) == 0.0)
      {
         ipn.array() = cs(3)*(cs(0)*ipn_2.array() + cs(1)*igrid.array()*ipn_1.array());
      }
      else
      {
         ipn.array() = cs(3)*(cs(0)*ipn_2.array() + (cs(1)*igrid.array() + cs(2))*ipn_1.array());
      }
   }

   void ThreeTermRecurrence::dPn(Eigen::Ref<internal::Matrix> idpn, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipn, const Eigen::Ref<const internal::Matrix>& ipn_1, const internal::Array& igrid, NormalizerNAB norm)
   {
      // Karniadakis and Sherwin page 586
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn, alpha, beta);
      idpn.array() =  ((cs(1) + cs(2)*igrid.array())*ipn.array()
         + cs(3)*ipn_1.array())
         / (cs(0)*(MHD_MP(1.0) - igrid.array().square()).array());
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<internal::Matrix> ip1, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ip0, const internal::Array& igrid, NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      if(cs(1) == 0.0)
      {
         ip1.array() = (cs(2)*cs(0))*igrid.array()*ip0.array();
      }
      else
      {
         ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1))*ip0.array();
      }
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<internal::Matrix> ip1, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1));
   }

   void ThreeTermRecurrence::P0(Eigen::Ref<internal::Matrix> ip0, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      ip0.setConstant(cs(0));
   }

   ThreeTermRecurrence::ThreeTermRecurrence()
   {
   }

   ThreeTermRecurrence::~ThreeTermRecurrence()
   {
   }

}
}
