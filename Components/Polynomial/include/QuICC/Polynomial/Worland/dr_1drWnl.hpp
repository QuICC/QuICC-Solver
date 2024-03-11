/**
 * @file dr_1drWnl.hpp
 * @brief Implementation of the D 1/r D r Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_DR_1DRWNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_DR_1DRWNL_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the D 1/r D r Worland polynomial
    */
   class dr_1drWnl: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         dr_1drWnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         dr_1drWnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta): WorlandBase(alpha, dBeta){};

         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void dr_1drWnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      using namespace Internal::Literals;
      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland D 1/r D r operator with l < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      Internal::MHDFloat a = this->alpha(l);
      Internal::MHDFloat b = this->beta(l);
      Internal::MHDFloat a1 = this->alpha(l) + 1.0_mp;
      Internal::MHDFloat b1 = this->beta(l) + 1.0_mp;
      Internal::MHDFloat a2 = this->alpha(l) + 2.0_mp;
      Internal::MHDFloat b2 = this->beta(l) + 2.0_mp;
      Internal::MHDFloat dl = Internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      Internal::Array ixgrid = 2.0_mp*igrid.array()*igrid.array() - 1.0_mp;

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      Internal::Matrix ipnab(gN,2);
      Internal::Matrix idpnab(gN,2);
      Internal::Matrix id2pnab(gN,2);

      // Compute P_0
      this->computeW0l(ipnab.col(0), l-2, a, b, igrid, WorlandBase::normWP0ab());
      ipnab.col(0) *= (dl - 1.0_mp)*(dl + 1.0_mp);
      if(scale.size() > 0)
      {
         ipnab.col(0).array() *= scale.array();
      }

      // Compute l P
      if(l == 1)
      {
         ipnab.col(1).setZero();
         evaluator(rOut, ipnab.col(1), 0);
      } else
      {
         evaluator(rOut, ipnab.col(0), 0);
      }

      if(nPoly > 1)
      {
         // Compute P_0
         if(l != 1)
         {
            ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid, WorlandBase::normWP1ab());
         }

         // Compute DP_1
         this->computeW0l(idpnab.col(0), l, a1, b1, igrid, WorlandBase::normWDP0ab());
         idpnab.col(0) *= 2.0_mp*(dl + 1.0_mp);
         if(scale.size() > 0)
         {
            idpnab.col(0).array() *= scale.array();
         }

         // Compute e P + 4r^2 DP
         if(l == 1)
         {
            evaluator(rOut, idpnab.col(0), 1);
         } else
         {
            evaluator(rOut, ipnab.col(1) + idpnab.col(0), 1);
         }
      }

      if(nPoly > 2)
      {
         if(l != 1)
         {
            // Increment P_n
            ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));
         }

         // Compute DP_1
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, WorlandBase::normWDP1ab());

         // Compute D2P_0
         this->computeW0l(id2pnab.col(0), l+2, a2, b2, igrid, WorlandBase::normWD2P0ab());
         if(scale.size() > 0)
         {
            id2pnab.col(0).array() *= scale.array();
         }

         // Compute e P + 2(x+1) DP
         if(l == 1)
         {
            evaluator(rOut, idpnab.col(1) + id2pnab.col(0), 2);
         } else
         {
            evaluator(rOut, ipnab.col(1) + idpnab.col(1) + id2pnab.col(0), 2);
         }
      }

      if(nPoly > 3)
      {
         if(l != 1)
         {
            // Increment P_3
            ThreeTermRecurrence::Pn(ipnab.col(0), 3, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));
         }

         // Compute DP_2
         ThreeTermRecurrence::Pn(idpnab.col(0), 2, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_1
         ThreeTermRecurrence::P1(id2pnab.col(1), a2, b2, id2pnab.col(0), ixgrid, WorlandBase::normWD2P1ab());

         // Compute e P + 2(x+1) DP
         if(l == 1)
         {
            evaluator(rOut, idpnab.col(1) + id2pnab.col(1), 3);
         } else
         {
            evaluator(rOut, ipnab.col(1) + idpnab.col(1) + id2pnab.col(1), 3);
         }
      }

      for(int i = 4; i < nPoly; ++i)
      {
         if(l != 1)
         {
            // Increment P_n
            ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));
         }

         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
         idpnab.col(0).swap(idpnab.col(1));

         // Increment D2P_n
         ThreeTermRecurrence::Pn(id2pnab.col(0), i-2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, WorlandBase::normWD2Pnab());
         id2pnab.col(0).swap(id2pnab.col(1));

         // Compute e P + 2(x+1) DP
         if(l == 1)
         {
            evaluator(rOut, idpnab.col(1) + id2pnab.col(1), i);
         } else
         {
            evaluator(rOut, ipnab.col(1) + idpnab.col(1) + id2pnab.col(1), i);
         }
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_DR_1DRWNL_HPP
