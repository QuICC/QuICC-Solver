/**
 * @file dr_1Wnl.hpp
 * @brief Implementation of d/dr[W_n^l / r] for Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_DR_1WNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_DR_1WNL_HPP

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
    * @brief Implementation of d/dr[W_n^l / r]
    *
    * The formula is:
    *   d/dr[W_n^l/r] = c_nl r^{l-2} [(l-1) P_n^{(-1/2, l-1/2)}
    *                                  + 2r^2(n+l) P_{n-1}^{(1/2, l+1/2)}]
    *
    * For l=1 the first term vanishes (l-1=0) so only the second term survives,
    * which is handled separately to avoid a r^{-1} seed that would be NaN at origin.
    * l=0 and l<0 are undefined and throw.
    */
   class dr_1Wnl: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         dr_1Wnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         dr_1Wnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta): WorlandBase(alpha, dBeta){};

         /**
          * @brief Compute d/dr[W_n^l / r] at grid points
          *
          * @tparam TEvaluator The evaluator allows to change behavior from computing Matrix operator, to On-the-fly transforms, etc
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

      protected:
         /**
          * @brief Special case for l = 1
          *
          * When l=1 the (l-1) prefactor of the first term is zero, leaving only
          * the derivative term: c_n1 * 2r(n+1) * P_{n-1}^{(1/2, 3/2)}.
          * This avoids seeding with r^{-1} which would be NaN at the origin.
          */
         template <typename T, typename TEvaluator> void computedr_1Wn1(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

      private:

   };

   template <typename T, typename TEvaluator> void dr_1Wnl::computedr_1Wn1(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      using namespace Internal::Literals;
      int gN = igrid.rows();

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not match grid size");
      }

      Internal::Matrix idiff(gN,2);

      // For l=1: (a+1, b+1) = (1/2, 3/2)
      Internal::MHDFloat a1 = this->alpha(1) + 1.0_mp;
      Internal::MHDFloat b1 = this->beta(1) + 1.0_mp;

      // Make X grid in [-1, 1]
      Internal::Array ixgrid = 2.0_mp*igrid.array()*igrid.array() - 1.0_mp;

      // n=0: d/dr[W_0^1/r] = 0  (P_{-1} = 0, first term has l-1=0)
      idiff.col(0).setZero();
      evaluator(rOut, idiff.col(0), 0);

      if(nPoly > 1)
      {
         // n=1: seed second-term recurrence at r^1 (= r^{l-2} * r^2 = r^{-1} * r^2)
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

   template <typename T, typename TEvaluator> void dr_1Wnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      using namespace Internal::Literals;
      if(l < 0)
      {
         throw std::logic_error("Tried to compute d/dr[W_n^l/r] with l < 0");
      }

      if(l == 0)
      {
         // First term needs to vanish. Artificially set it to zero.
         //tbdone 
         //this->computedr_1Wn0(rOut, nPoly, igrid, scale, evaluator);
      }
      else if(l == 1)
      {
         // First term vanishes (l-1=0); only the derivative term survives.
         this->computedr_1Wn1(rOut, nPoly, igrid, scale, evaluator);
      } else
      {
         // General case l >= 2
         int gN = igrid.rows();

         if (nPoly < 1)
         {
            throw std::logic_error("Operator matrix should have at least 1 column");
         }

         if (gN != igrid.size())
         {
            throw std::logic_error("Operator matrix does not match grid size");
         }

         Internal::MHDFloat a = this->alpha(l);
         Internal::MHDFloat b = this->beta(l);
         Internal::MHDFloat a1 = this->alpha(l) + 1.0_mp;
         Internal::MHDFloat b1 = this->beta(l) + 1.0_mp;
         // Prefactor of first term is (l-1), not l
         Internal::MHDFloat dl = Internal::MHDFloat(l - 1);

         // Make X grid in [-1, 1]
         Internal::Array ixgrid = 2.0_mp*igrid.array()*igrid.array() - 1.0_mp;

         // Storage for P_n^{(alpha,beta)} and dP_n^{(alpha+1,beta+1)}
         Internal::Matrix ipnab(gN,2);
         Internal::Matrix idpnab(gN,2);

         // Seed first term at r^{l-2}  (output lives at Worland level l-2)
         this->computeW0l(ipnab.col(0), l-2, a, b, igrid, WorlandBase::normWP0ab());
         ipnab.col(0) *= dl;
         if(scale.size() > 0)
         {
            ipnab.col(0).array() *= scale.segment(0,gN).array();
         }

         // Second term: P_{-1} = 0 at n=0
         idpnab.col(0).setZero();

         // Compute (l-1) r^{l-2} P_0
         evaluator(rOut, ipnab.col(0), 0);

         if(nPoly > 1)
         {
            // Advance first-term recurrence to P_1
            ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid, WorlandBase::normWP1ab());

            // Seed second term at r^l  (= r^{l-2} * r^2)
            this->computeW0l(idpnab.col(0), l, a1, b1, igrid, WorlandBase::normWDP0ab());
            if(scale.size() > 0)
            {
               idpnab.col(0).array() *= scale.segment(0,gN).array();
            }

            // (l-1) r^{l-2} P_1 + r^l c_DP P_0^{(a1,b1)}
            evaluator(rOut, ipnab.col(1) + idpnab.col(0), 1);
         }

         if(nPoly > 2)
         {
            // Advance first-term recurrence to P_2
            ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));

            // Advance second-term recurrence to P_1^{(a1,b1)}
            ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, WorlandBase::normWDP1ab());

            evaluator(rOut, ipnab.col(1) + idpnab.col(1), 2);
         }

         for(int i = 3; i < nPoly; ++i)
         {
            // Advance first-term recurrence
            ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0), ixgrid, WorlandBase::normWPnab());
            ipnab.col(0).swap(ipnab.col(1));

            // Advance second-term recurrence
            ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
            idpnab.col(0).swap(idpnab.col(1));

            // (l-1) r^{l-2} P_n + r^l c_DP P_{n-1}^{(a1,b1)}
            evaluator(rOut, ipnab.col(1) + idpnab.col(1), i);
         }
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_DR_1WNL_HPP
