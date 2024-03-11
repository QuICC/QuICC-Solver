/**
 * @file r_1Wnl.hpp
 * @brief Implementation of the 1/r Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_R_1WNLEXPLICIT_HPP
#define QUICC_POLYNOMIAL_WORLAND_R_1WNLEXPLICIT_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /// @brief Generic implementation
   /// @tparam
   template <class >
   class r_1Wnl;

   /**
    * @brief Implementation of the Worland polynomial with explcit grid division
    */
   template <>
   class r_1Wnl<explicit_t>: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         r_1Wnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          */
         r_1Wnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta): WorlandBase(alpha, dBeta){};

         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEval> inline void r_1Wnl<explicit_t>::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEval evaluator)
   {
      Internal::Matrix tOp(igrid.size(), nPoly);
      Polynomial::Worland::Wnl wnl(this->alpha(l), this->dBeta());
      wnl.compute<Internal::MHDFloat>(tOp, nPoly, l, igrid, scale, evaluator);

      rOut = ((tOp.transpose()*igrid.array().pow(-1).matrix().asDiagonal()).transpose()).cast<T>();
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_R_1WNLEXPLICIT_HPP
