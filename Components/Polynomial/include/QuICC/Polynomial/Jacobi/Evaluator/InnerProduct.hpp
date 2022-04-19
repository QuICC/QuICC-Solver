/** 
 * @file InnerProduct.hpp
 * @brief Evaluator to compute inner product with polynmial
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_EVALUATOR_INNERPRODUCT_HPP
#define QUICC_POLYNOMIAL_JACOBI_EVALUATOR_INNERPRODUCT_HPP

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

namespace Jacobi {

namespace Evaluator {

   /**
    * @brief Evaluator to compute inner product with polynomial
    */ 
   template <typename T> class InnerProduct
   {
      public:
         /**
          * @brief Constructor
          */
         InnerProduct(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& in);

         /**
          * @brief Apply evaluator
          */
         void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i);

      protected:

      private:
         /**
          * @brief Reference to input expression
          */
         const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mIn;
   };

   template <typename T> InnerProduct<T>::InnerProduct(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& in)
      : mIn(in)
   {
   }

   template <typename T> inline void InnerProduct<T>::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i)
   {
      rOut.row(i) = Precision::cast(ipolycol).transpose()*this->mIn;
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_JACOBI_EVALUATOR_INNERPRODUCT_HPP
