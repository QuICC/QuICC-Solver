/** 
 * @file Set.hpp
 * @brief Evaluator to create matrix operator of polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_EVALUATOR_SET_HPP
#define QUICC_POLYNOMIAL_WORLAND_EVALUATOR_SET_HPP

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

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Evaluator {

   /**
    * @brief Evaluator to create matrix operator of polynomial
    */ 
   class Set
   {
      public:
         /**
          * @brief Set output to current column, casting to double precision if required
          */
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i);
   };

   template <typename T> inline void Set::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i)
   {
      rOut.col(i) = Precision::cast(ipolycol);
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_EVALUATOR_SET_HPP
