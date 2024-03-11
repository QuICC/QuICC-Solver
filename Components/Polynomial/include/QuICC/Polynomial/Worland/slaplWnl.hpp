/**
 * @file slaplWnl.hpp
 * @brief Implementation of the spherical laplacian of Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_SLAPLWNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_SLAPLWNL_HPP

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
    * @brief Implementation of the spherical laplacian of Worland polynomial
    */
   class slaplWnl: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         slaplWnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         slaplWnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta): WorlandBase(alpha, dBeta){};

         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEvaluator> void slaplWnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      using namespace Internal::Literals;
      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland spherical laplacian with l < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      Internal::MHDFloat a1 = this->alpha(l) + 1.0_mp;
      Internal::MHDFloat b1 = this->beta(l) + 1.0_mp;
      Internal::MHDFloat a2 = this->alpha(l) + 2.0_mp;
      Internal::MHDFloat b2 = this->beta(l) + 2.0_mp;
      Internal::MHDFloat dl = Internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      Internal::Array ixgrid = 2.0_mp*igrid.array()*igrid.array() - 1.0_mp;

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      Internal::Matrix idpnab(gN,2);
      Internal::Matrix id2pnab(gN,2);

      // Compute spherical laplacian P_0
      idpnab.col(0).setZero();
      evaluator(rOut, idpnab.col(0), 0);

      if(nPoly > 1)
      {
         // Compute DP_1
         this->computeW0l(idpnab.col(0), l, a1, b1, igrid, WorlandBase::normWDP0ab());
         idpnab.col(0) *= (2.0_mp*dl + 3.0_mp);
         if(scale.size() > 0)
         {
            idpnab.col(0).array() *= scale.segment(0,gN).array();
         }

         // Compute spherical laplacian P_1
         evaluator(rOut, idpnab.col(0), 1);
      }

      if(nPoly > 2)
      {
         // Increment DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, WorlandBase::normWDP1ab());

         // Compute D2P_2
         this->computeW0l(id2pnab.col(0), l+2, a2, b2, igrid, WorlandBase::normWD2P0ab());
         if(scale.size() > 0)
         {
            id2pnab.col(0).array() *= scale.segment(0,gN).array();
         }

         // Compute e P + 2(x+1) DP
         evaluator(rOut, id2pnab.col(0) + idpnab.col(1), 2);
      }

      if(nPoly > 3)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), 2, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_2
         ThreeTermRecurrence::P1(id2pnab.col(1), a2, b2, id2pnab.col(0), ixgrid, WorlandBase::normWD2P1ab());

         // Compute e P + 2(x+1) DP
         evaluator(rOut, id2pnab.col(1) + idpnab.col(1), 3);
      }

      for(int i = 4; i < nPoly; ++i)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
         idpnab.col(0).swap(idpnab.col(1));

         // Increment D2P_n
         ThreeTermRecurrence::Pn(id2pnab.col(0), i-2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, WorlandBase::normWD2Pnab());
         id2pnab.col(0).swap(id2pnab.col(1));

         // Compute e P + 2(x+1) DP
         evaluator(rOut, id2pnab.col(1) + idpnab.col(1), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_SLAPLWNL_HPP
