/**
 * @file Set.hpp
 * @brief Evaluator to compute matrix operator of polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_SET_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_SET_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Casts.hpp"

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
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i);

         /**
          * @brief Apply evaluator assuming even/odd symmetry
          */
         template <typename T> void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i, const bool isEven);
   };

   template <typename T> inline void Set::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i)
   {
      rOut.col(i) = Internal::cast(ipolycol);
   }

   template <typename T> inline void Set::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i, const bool isEven)
   {
      int gN = ipolycol.rows();

      rOut.col(i).topRows(gN) = Internal::cast(ipolycol);

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
