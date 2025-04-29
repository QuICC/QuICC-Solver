/**
 * @file OuterProduct.hpp
 * @brief Evaluator to compute outer product with polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_EVALUATOR_OUTERPRODUCT_HPP
#define QUICC_POLYNOMIAL_WORLAND_EVALUATOR_OUTERPRODUCT_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Casts.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Evaluator {

   /**
    * @brief Evaluator to compute outer product with polynomial
    */
   template <typename T> class OuterProduct
   {
      public:
         /**
          * @brief Constructor
          */
         OuterProduct(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& in);

         /**
          * @brief Apply evaluator
          */
         void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i);

      protected:

      private:
         /**
          * @brief Reference to input expression
          */
         const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mIn;
   };

   template <typename T> OuterProduct<T>::OuterProduct(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& in)
      : mIn(in)
   {
   }

   template <typename T> inline void OuterProduct<T>::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i)
   {
      if(i == 0)
      {
         rOut = Internal::cast(ipolycol)*this->mIn.row(i);
      } else
      {
         rOut += Internal::cast(ipolycol)*this->mIn.row(i);
      }
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_EVALUATOR_OUTERPRODUCT_HPP
