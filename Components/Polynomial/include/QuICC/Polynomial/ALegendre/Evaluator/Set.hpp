/** 
 * @file Set.hpp
 * @brief Evaluator to compute matrix operator of polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_SET_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_SET_HPP

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
    * @brief Evaluator to compute matrix operator of polynomial
    */ 
   class Set
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

   template <typename T> inline void Set::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i)
   {  
      rOut.col(i) = Precision::cast(ipolycol);
   }

   template <typename T> inline void Set::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i, const bool isEven)
   {  
      int gN = ipolycol.rows();

      rOut.col(i).topRows(gN) = Precision::cast(ipolycol);

      if(isEven)
      {
         rOut.col(i).bottomRows(gN) = rOut.col(i).topRows(gN).reverse();
      } else
      {
         rOut.col(i).bottomRows(gN) = -rOut.col(i).topRows(gN).reverse();
      }
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_SET_HPP
