/**
 * @file claplhWnl.hpp
 * @brief Implementation of the cylindrical laplacian Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_CLAPLHWNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_CLAPLHWNL_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Precision.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the cylindrical laplacian Worland polynomial
    */
   class claplhWnl: public WorlandBase
   {
      public:
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEvaluator> void claplhWnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland cylindrical laplacian with l < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a1 = this->alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = this->beta(l) + MHD_MP(1.0);
      internal::MHDFloat a2 = this->alpha(l) + MHD_MP(2.0);
      internal::MHDFloat b2 = this->beta(l) + MHD_MP(2.0);
      internal::MHDFloat dl = internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix idpnab(gN,2);
      internal::Matrix id2pnab(gN,2);

      // Compute spherical laplacian P_0
      idpnab.col(0).setZero();
      evaluator(rOut, idpnab.col(0), 0);

      if(nPoly > 1)
      {
         // Compute DP_1
         this->computeW0l(idpnab.col(0), l, a1, b1, igrid, WorlandBase::normWDP0ab());
         idpnab.col(0) *= MHD_MP(2.0)*(dl + MHD_MP(1.0));
         if(scale.size() > 0)
         {
            idpnab.col(0).array() *= scale.array();
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
            id2pnab.col(0).array() *= scale.array();
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

#endif // QUICC_POLYNOMIAL_WORLAND_CLAPLHWNL_HPP
