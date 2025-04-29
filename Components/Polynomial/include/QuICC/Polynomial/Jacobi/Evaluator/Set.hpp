/**
 * @file Set.hpp
 * @brief Evaluator to create matrix operator of polynomial
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_EVALUATOR_SET_HPP
#define QUICC_POLYNOMIAL_JACOBI_EVALUATOR_SET_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Casts.hpp"

namespace QuICC {

namespace Polynomial {

namespace Jacobi {

namespace Evaluator {

   /**
    * @brief Evaluator to create matrix operator of polynomial
    */
   class Set
   {
      public:
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i);
   };

   template <typename T> inline void Set::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i)
   {
      rOut.col(i) = Internal::cast(ipolycol);
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_JACOBI_EVALUATOR_SET_HPP
