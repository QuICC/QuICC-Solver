/**
 * @file Wnl.hpp
 * @brief Implementation of the Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_WNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_WNL_HPP

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
    * @brief Implementation of the Worland polynomial
    */
   class Wnl: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         Wnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param lShift  Shift in harmonic degree
          */
         Wnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta, const int lShift = 0): WorlandBase(alpha, dBeta), mLShift(lShift){};

         /**
          * @brief Compute worland polynomial
          *
          * @tparam TEvaluator The evaluator allows to change behavior from computing Matric operator, to On-the-fly transforms, etc
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

      private:
         /**
          * @brief Shift of the harmonic degree when computing the polyniomial
          */
         int mLShift = 0;
   };

   template <typename T, typename TEval> inline void Wnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int lIn, const Internal::Array& igrid, const Internal::Array& scale, TEval evaluator)
   {
      using namespace Internal::Literals;
      const auto l = lIn + this->mLShift;

      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland polynomial W_n^l with l < 0");
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

      Internal::MHDFloat a = this->alpha(l);
      Internal::MHDFloat b = this->beta(l);

      this->computeW0l(ipoly.col(0), l, a, b, igrid, WorlandBase::normWP0ab());
      if(scale.size() > 0)
      {
         ipoly.col(0).array() *= scale.array();
      }
      evaluator(rOut, ipoly.col(0), 0);

      // Make X grid in [-1, 1]
      Internal::Array ixgrid = 2.0_mp*igrid.array()*igrid.array() - 1.0_mp;

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), a, b, ipoly.col(0), ixgrid, WorlandBase::normWP1ab());
         evaluator(rOut, ipoly.col(1), 1);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(0), i, a, b, ipoly.col(1), ipoly.col(0), ixgrid, WorlandBase::normWPnab());
         ipoly.col(0).swap(ipoly.col(1));
         evaluator(rOut, ipoly.col(1), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_WNL_HPP
