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

   void ThreeTermRecurrence::Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerNC norm)
   {
      Internal::MHDFloat dn = Internal::MHDFloat(n);
      Internal::Array cs = norm(dn);

      ipn.array() = cs(3)*(cs(0)*ipn_2.array() + (cs(1)*igrid.array() + cs(2))*ipn_1.array());
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<Internal::Matrix> ip1, const Eigen::Ref<const Internal::Matrix>& ip0, const Internal::Array& igrid, NormalizerC norm)
   {
      Internal::Array cs = norm();

      ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1))*ip0.array();
   }

   void ThreeTermRecurrence::P0(Eigen::Ref<Internal::Matrix> ip0, const Internal::Array& igrid, NormalizerC norm)
   {
      Internal::Array cs = norm();

      ip0.setConstant(cs(0));
   }

   void ThreeTermRecurrence::Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Internal::MHDFloat alpha, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerNA norm)
   {
      Internal::MHDFloat dn = Internal::MHDFloat(n);
      Internal::Array cs = norm(dn, alpha);

      ipn.array() = cs(3)*(cs(0)*ipn_2.array() + (cs(1)*igrid.array() + cs(2))*ipn_1.array());
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::MHDFloat alpha, const Eigen::Ref<const Internal::Matrix>& ip0, const Internal::Array& igrid, NormalizerA norm)
   {
      Internal::Array cs = norm(alpha);

      ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1))*ip0.array();
   }

   void ThreeTermRecurrence::P0(Eigen::Ref<Internal::Matrix> ip0, const Internal::MHDFloat alpha, const Internal::Array& igrid, NormalizerA norm)
   {
      Internal::Array cs = norm(alpha);

      ip0.setConstant(cs(0));
   }

   void ThreeTermRecurrence::Pn(Eigen::Ref<Internal::Matrix> ipn, const int n, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Eigen::Ref<const Internal::Matrix>& ipn_2, const Internal::Array& igrid, NormalizerNAB norm)
   {
      Internal::MHDFloat dn = Internal::MHDFloat(n);
      Internal::Array cs = norm(dn, alpha, beta);

      if(cs(2) == 0.0)
      {
         ipn.array() = cs(3)*(cs(0)*ipn_2.array() + cs(1)*igrid.array()*ipn_1.array());
      }
      else
      {
         ipn.array() = cs(3)*(cs(0)*ipn_2.array() + (cs(1)*igrid.array() + cs(2))*ipn_1.array());
      }
   }

   void ThreeTermRecurrence::dPn(Eigen::Ref<Internal::Matrix> idpn, const int n, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Eigen::Ref<const Internal::Matrix>& ipn, const Eigen::Ref<const Internal::Matrix>& ipn_1, const Internal::Array& igrid, NormalizerNAB norm)
   {
      // Karniadakis and Sherwin page 586
      Internal::MHDFloat dn = Internal::MHDFloat(n);
      Internal::Array cs = norm(dn, alpha, beta);
      idpn.array() =  ((cs(1) + cs(2)*igrid.array())*ipn.array()
         + cs(3)*ipn_1.array())
         / (cs(0)*(MHD_MP(1.0) - igrid.array().square()).array());
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Eigen::Ref<const Internal::Matrix>& ip0, const Internal::Array& igrid, NormalizerAB norm)
   {
      Internal::Array cs = norm(alpha, beta);

      if(cs(1) == 0.0)
      {
         ip1.array() = (cs(2)*cs(0))*igrid.array()*ip0.array();
      }
      else
      {
         ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1))*ip0.array();
      }
   }

   void ThreeTermRecurrence::P1(Eigen::Ref<Internal::Matrix> ip1, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, NormalizerAB norm)
   {
      Internal::Array cs = norm(alpha, beta);

      ip1.array() = cs(2)*(cs(0)*igrid.array() + cs(1));
   }

   void ThreeTermRecurrence::P0(Eigen::Ref<Internal::Matrix> ip0, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, NormalizerAB norm)
   {
      Internal::Array cs = norm(alpha, beta);

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
