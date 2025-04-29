/**
 * @file r_1drWnlExplicit.hpp
 * @brief Implementation of the 1/r D r Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_R_1DRWNLEXPLICIT_HPP
#define QUICC_POLYNOMIAL_WORLAND_R_1DRWNLEXPLICIT_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /// @brief Generic implementation
   /// @tparam
   template <class >
   class r_1drWnl;

   /**
    * @brief Implementation of the Worland polynomial with explicit grid division
    */
   template <>
   class r_1drWnl<explicit_t>: public WorlandBase
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

   template <typename T, typename TEval> inline void r_1drWnl<explicit_t>::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEval evaluator)
   {
      // assert(false && "this is not working properly");
      if (scale.size() == 0)
      {
         throw std::logic_error("this polynomial implementation needs to know the weights for the internal projection");
      }

      if(igrid.size() < (2*nPoly + l + 2)/2)
      {
         throw std::logic_error("Grid size does not allow exact integration for internal step");
      }

      Internal::Matrix tOp(igrid.size(), nPoly);

      // Extend intermediate truncation by one due to multiplication by r
      Internal::Matrix opA(igrid.size(), nPoly+1);
      Polynomial::Worland::Wnl wnl(this->alpha(l), this->dBeta());
      wnl.compute<Internal::MHDFloat>(opA, nPoly+1, l, igrid, scale.array()*igrid.array(), evaluator);

      Internal::Matrix opB(igrid.size(), nPoly+1);
      Polynomial::Worland::dWnl dWnl(this->alpha(l), this->dBeta());
      dWnl.compute<Internal::MHDFloat>(opB, nPoly+1, l, igrid, Internal::Array(), evaluator);

      Internal::Matrix opC(igrid.size(), nPoly);
      wnl.compute<Internal::MHDFloat>(opC, nPoly, l, igrid, scale.array()*igrid.array().pow(-1), evaluator);

      tOp = (opC.transpose()*opB*opA.transpose()).transpose();

      rOut = tOp.cast<T>();
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_R_1DRWNLEXPLICIT_HPP
