/** 
 * @file dWnl.hpp
 * @brief Implementation of the first derivative of Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_DWNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_DWNL_HPP

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
#include "QuICC/Precision.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the first derivative of Worland polynomial
    */ 
   class dWnl: public WorlandBase
   {
      public:
         /**
          * @brief Compute worland polynomial
          *
          * @tparam TEValuator The evaluator allows to change behavior from computing Matric operator, to On-the-fly transforms, etc
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);

      protected:
         /**
          * @brief Special case for l = 0
          */
         template <typename T, typename TEvaluator> void computedWn0(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);

      private:

   };

   template <typename T, typename TEvaluator> void dWnl::computedWn0(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      internal::Matrix idiff(gN,2);

      internal::MHDFloat a1 = this->alpha(0) + MHD_MP(1.0);
      internal::MHDFloat b1 = this->beta(0) + MHD_MP(1.0);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      idiff.col(0).setZero();
      evaluator(rOut, idiff.col(0), 0);

      if(nPoly > 1)
      {
         this->computeW0l(idiff.col(1), 1, a1, b1, igrid, WorlandBase::normWDP0ab());
         if(scale.size() > 0)
         {
            idiff.col(1).array() *= scale.segment(0,gN).array();
         }
         evaluator(rOut, idiff.col(1), 1);
      }

      if(nPoly > 2)
      {
         ThreeTermRecurrence::P1(idiff.col(0), a1, b1, idiff.col(1), ixgrid, WorlandBase::normWDP1ab());
         idiff.col(0).swap(idiff.col(1));
         evaluator(rOut, idiff.col(1), 2);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(0), i-1, a1, b1, idiff.col(1), idiff.col(0), ixgrid, WorlandBase::normWDPnab());
         idiff.col(0).swap(idiff.col(1));
         evaluator(rOut, idiff.col(1), i);
      }
   }

   template <typename T, typename TEvaluator> void dWnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
   {
      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland derivative with l < 0");
      }

      if(l == 0)
      {
         this->computedWn0(rOut, nPoly, igrid, scale, evaluator);
      } else
      {
         int gN = igrid.rows();

         if (nPoly < 1)
         {
            throw std::logic_error("Operator matrix should have at least 1 column");
         }

         if (gN != igrid.size())
         {
            throw std::logic_error("Operator matrix does not mach grid size");
         }

         internal::MHDFloat a = this->alpha(l);
         internal::MHDFloat b = this->beta(l);
         internal::MHDFloat a1 = this->alpha(l) + MHD_MP(1.0);
         internal::MHDFloat b1 = this->beta(l) + MHD_MP(1.0);
         internal::MHDFloat dl = internal::MHDFloat(l);

         // Make X grid in [-1, 1]
         internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

         // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
         internal::Matrix ipnab(gN,2);
         internal::Matrix idpnab(gN,2);

         // Compute P_0
         this->computeW0l(ipnab.col(0), l-1, a, b, igrid, WorlandBase::normWP0ab());
         ipnab.col(0) *= dl;
         if(scale.size() > 0)
         {
            ipnab.col(0).array() *= scale.segment(0,gN).array();
         }

         // Compute DP_0
         idpnab.col(0).setZero();

         // Compute l P
         evaluator(rOut, ipnab.col(0), 0);

         if(nPoly > 1)
         {
            // Compute P_0
            ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid, WorlandBase::normWP1ab());

            // Compute DP_1
            this->computeW0l(idpnab.col(0), l+1, a1, b1, igrid, WorlandBase::normWDP0ab());
            if(scale.size() > 0)
            {
               idpnab.col(0).array() *= scale.segment(0,gN).array();
            }

            // Compute e P + 4r^2 DP
            evaluator(rOut, ipnab.col(1) + idpnab.col(0), 1);
         }

         if(nPoly > 2)
         {
            // Increment P_n
            ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));

            // Compute DP_2
            ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, WorlandBase::normWDP1ab());

            // Compute e P + 2(x+1) DP
            evaluator(rOut, ipnab.col(1) + idpnab.col(1), 2);
         }

         for(int i = 3; i < nPoly; ++i)
         {
            // Increment P_n
            ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));

            // Increment DP_n
            ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
            idpnab.col(0).swap(idpnab.col(1));

            // Compute e P + 2(x+1) DP
            evaluator(rOut, ipnab.col(1) + idpnab.col(1), i);
         }
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_DWNL_HPP
