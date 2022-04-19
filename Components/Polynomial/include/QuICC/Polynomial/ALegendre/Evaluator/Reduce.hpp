/** 
 * @file Reduce.hpp
 * @brief Evaluator to compute column sum reduction of polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_REDUCE_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_REDUCE_HPP

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

namespace ALegendre {

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
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i);

         /**
          * @brief Apply evaluator assuming even/odd symmetry
          */
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i, const bool isEven);
   };

   template <typename T> void Reduce::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i)
   {
      rOut(i,0) = Precision::cast(ipolycol.sum());
   }

   template <typename T> void Reduce::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i, const bool isEven)
   {
      if(isEven)
      {
         rOut(i,0) = 2.0*Precision::cast(ipolycol.sum());
      } else
      {
         rOut(i,0) = 0.0;
      }
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_REDUCE_HPP
