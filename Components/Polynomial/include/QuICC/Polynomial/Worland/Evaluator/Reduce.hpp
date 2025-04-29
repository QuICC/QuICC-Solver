/**
 * @file Reduce.hpp
 * @brief Evaluator to compute column sum reduction of polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_EVALUATOR_REDUCE_HPP
#define QUICC_POLYNOMIAL_WORLAND_EVALUATOR_REDUCE_HPP

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
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Evaluator {

   /**
    * @brief Evaluator to compute column sum reduction of polynomial
    */
   class Reduce
   {
      public:
         /**
          * @brief Apply evaluator
          */
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i);
   };

   template <typename T> inline void Reduce::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i)
   {
      rOut(i,0) = Internal::cast(ipolycol.sum());
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_EVALUATOR_REDUCE_HPP
