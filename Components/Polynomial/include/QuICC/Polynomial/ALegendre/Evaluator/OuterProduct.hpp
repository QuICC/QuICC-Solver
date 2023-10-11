/**
 * @file OuterProduct.hpp
 * @brief Evaluator to compute outer product with polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_OUTERPRODUCT_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_OUTERPRODUCT_HPP

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

namespace ALegendre {

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
          * @brief Prepare grid for polynomial calculation
          */
         Internal::Array prepareGrid(const Internal::Array& fullGrid) const;

         /**
          * @brief Apply evaluator
          */
         void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i);

         /**
          * @brief Apply operation assuming even or odd symmetry
          */
         void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i, const bool isEven);

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

   template <typename T> inline void OuterProduct<T>::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const Internal::Matrix>& ipolycol, const int i, const bool isEven)
   {
      int gN = ipolycol.rows();

      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp = Internal::cast(ipolycol)*this->mIn.row(i);
      if(i == 0)
      {
         if(isEven)
         {
            rOut.bottomRows(gN) = tmp.colwise().reverse();
         } else
         {
            rOut.bottomRows(gN) = -tmp.colwise().reverse();
         }
         rOut.topRows(gN) = tmp;
      } else
      {
         if(rOut.rows()%2 == 1)
         {
            tmp.bottomRows(1) *= 0.5;
         }
         rOut.topRows(gN) += tmp;
         if(isEven)
         {
            rOut.bottomRows(gN) += tmp.colwise().reverse();
         } else
         {
            rOut.bottomRows(gN) -= tmp.colwise().reverse();
         }
      }
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_EVALUATOR_OUTERPRODUCT_HPP
