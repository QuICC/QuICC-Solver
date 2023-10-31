/**
 * @file rWnl.hpp
 * @brief Implementation of the Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_RWNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_RWNL_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the Worland polynomial
    */
   class rWnl: public WorlandBase
   {
      public:
         /**
          * @brief Default constructor
          */
         rWnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          */
         rWnl(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta);

         /**
          * @brief Compute worland polynomial
          *
          * @tparam TEvaluator The evaluator allows to change behavior from computing Matric operator, to On-the-fly transforms, etc
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEval> inline void rWnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale, TEval evaluator)
   {
        Polynomial::Worland::Wnl wnl;

        wnl.compute<T>(rOut, nPoly, l, igrid, (scale.array()*igrid.array()).matrix(), evaluator);
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_RWNL_HPP
