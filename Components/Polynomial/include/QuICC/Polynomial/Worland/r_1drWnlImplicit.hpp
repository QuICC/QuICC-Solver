/**
 * @file r_1drWnlImplicit.hpp
 * @brief Implementation of the 1/r D r Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_R_1DRWNLIMPLICIT_HPP
#define QUICC_POLYNOMIAL_WORLAND_R_1DRWNLIMPLICIT_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /// @brief Generic implementation
   /// @tparam
   template <class >
   class r_1drWnl;

   /**
    * @brief Implementation of the Worland polynomial with implicit grid division
    */
   template <>
   class r_1drWnl<implicit_t>: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         r_1drWnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         r_1drWnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta): WorlandBase(alpha, dBeta){};

         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEval> inline void r_1drWnl<implicit_t>::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEval evaluator)
   {
      Polynomial::Worland::r_1drWnl<QuICC::Polynomial::Worland::recurrence_t> r_1drWnl(this->alpha(l), this->dBeta());

      if (l == 0)
      {
         // fallback on explicit implementation
         Polynomial::Worland::r_1drWnl<QuICC::Polynomial::Worland::explicit_t> r_1drWnlExplicit(this->alpha(l), this->dBeta());
         r_1drWnlExplicit.compute<T>(rOut, nPoly, l, igrid, scale, evaluator);
      }
      else
      {
         if (scale.size() == 0)
         {
            throw std::logic_error("this polynomial implementation needs to know the weights for the internal projection");
         }

         if(igrid.size() < (2*nPoly + l + 1)/2)
         {
            throw std::logic_error("Grid size does not allow exact integration for internal step");
         }

         // Operates on polynomials with l = l-1
         int lm = std::abs(l-1);
         int n_in = nPoly + 1;

         Internal::Matrix opA(igrid.size(), n_in);
         Polynomial::Worland::Wnl wnl(this->alpha(l), this->dBeta());
         wnl.compute<Internal::MHDFloat>(opA, n_in, lm, igrid, scale, evaluator);

         Internal::Matrix opB(igrid.size(), n_in);
         r_1drWnl.compute<Internal::MHDFloat>(opB, n_in, lm, igrid, Internal::Array(), evaluator);

         Internal::Matrix opC(igrid.size(), nPoly);
         Polynomial::Worland::Wnl wnlB(this->alpha(l), this->dBeta());
         wnlB.compute<Internal::MHDFloat>(opC, nPoly, l, igrid, scale, evaluator);

         rOut = ((opC.transpose()*opB*opA.transpose()).transpose()).cast<T>();
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_R_1DRWNLIMPLICIT_HPP
