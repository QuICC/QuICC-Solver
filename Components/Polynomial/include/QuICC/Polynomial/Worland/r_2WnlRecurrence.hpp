/**
 * @file r_2Wnl.hpp
 * @brief Implementation of the 1/r^2 Worland polynomial
 * @brief sets the leading order term in the l=1 case to zero
 */

 // TODO:
 // Improve compute function (see below)
 
#ifndef QUICC_POLYNOMIAL_WORLAND_R_2WNLRECURRENCE_HPP
#define QUICC_POLYNOMIAL_WORLAND_R_2WNLRECURRENCE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /// @brief Generic implementation
   /// @tparam
   template <class >
   class r_2Wnl;

   /**
    * @brief Implementation of the Worland polynomial with recurrence relation
    */
   template <>
   class r_2Wnl<recurrence_t>: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         r_2Wnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         r_2Wnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta): WorlandBase(alpha, dBeta){};

         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void r_2Wnl<recurrence_t>::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {

      // TODO:
      // For l=1 this code removes constant Jacobi term to eliminate the 1/r singularity.
      // Mathematically correct but computes (P(x)-P(-1))/r explicitly,
      // relying on evaluations of 1/r, potentailly bad near the origin.
      // A future improvement is to evaluate
      //     (P(x)-P(-1))/(x+1)
      // directly via a dedicated recurrence.

      using namespace Internal::Literals;
      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland 1/r operator with l < 0");
      }

      if(nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if(gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      Internal::Matrix ipoly(gN,2);

      // for l=1 we need to remove the leading order term (singular)
      Internal::Matrix ipoly0(gN,2);
      ipoly0.setZero();

      Internal::Matrix ipolyCorr(gN,2);
      ipolyCorr.setZero();

      Internal::MHDFloat a = this->alpha(l);
      Internal::MHDFloat b = this->beta(l);

      // Initialize
      this->computeW0l(ipoly.col(0), l-2, a, b, igrid, WorlandBase::normWP0ab());
      if(l==1)
      {
         // P_n^(a,b) in x = -1
         this->computeW0l(ipoly0.col(0), 0, a, b, 0*igrid, WorlandBase::normWP0ab());
      }

      if(scale.size() > 0)
      {
         ipoly.col(0).array() *= scale.array();
         ipoly0.col(0).array() *= scale.array();
      }
      // for l=1, remove P_n^(a,b)(x=-1)/r (which is exactly the singular term)
      ipolyCorr.col(0) = ipoly.col(0)-(igrid.array().pow(-1) * ipoly0.col(0).array()).matrix();
      evaluator(rOut, ipolyCorr.col(0), 0);

      // Make X grid in [-1, 1]
      Internal::Array ixgrid = 2.0_mp*igrid.array()*igrid.array() - 1.0_mp;
      // Make X "grid" of x=-1
      Internal::Array ixgrid0 = 0*igrid.array() - 1.0_mp;

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), a, b, ipoly.col(0), ixgrid, WorlandBase::normWP1ab());
         if(l==1)
         {
            ThreeTermRecurrence::P1(ipoly0.col(1), a, b, ipoly0.col(0), ixgrid0, WorlandBase::normWP1ab());
         }
         ipolyCorr.col(1) = ipoly.col(1)-(igrid.array().pow(-1) * ipoly0.col(1).array()).matrix();
         evaluator(rOut, ipolyCorr.col(1), 1);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(0), i, a, b, ipoly.col(1), ipoly.col(0), ixgrid, WorlandBase::normWPnab());
         ipoly.col(0).swap(ipoly.col(1));

         if(l==1)
         {
            ThreeTermRecurrence::Pn(ipoly0.col(0), i, a, b, ipoly0.col(1), ipoly0.col(0), ixgrid0, WorlandBase::normWPnab());
            ipoly0.col(0).swap(ipoly0.col(1));
         }

         ipolyCorr.col(1) = ipoly.col(1)-(igrid.array().pow(-1) * ipoly0.col(1).array()).matrix();
         evaluator(rOut, ipolyCorr.col(1), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_R_2WNLRECURRENCE_HPP
