/** 
 * @file InnerProduct.hpp
 * @brief Implementation of a simple on-the-fly function to set matrix column
 */

#ifndef QUICC_POLYNOMIAL_LEGENDRE_EVALUATOR_INNERPRODUCT_HPP
#define QUICC_POLYNOMIAL_LEGENDRE_EVALUATOR_INNERPRODUCT_HPP

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
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Polynomial {

namespace Legendre {

namespace Evaluator {

   /**
    * @brief Implementation of the Legendre polynomial
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
         inline void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i);

         /**
          * @brief Apply evaluator assuming even or odd symetry
          */
         inline void operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i, const bool isEven);

      protected:

      private:
         /**
          * @brief Reference to input expression
          */
         const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& mIn;

         /**
          * @brief Input storage with symmetry
          */
         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mEven;
         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mOdd;
   };

   template <typename T> InnerProduct<T>::InnerProduct(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >& in)
      : mIn(in)
   {
      int gN = in.rows()/2 + in.rows()%2;
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp = in.bottomRows(gN).colwise().reverse();
      this->mEven = in.topRows(gN) + tmp;
      this->mOdd = in.topRows(gN) - tmp;

      // Correct for odd number of grid points
      if(in.rows()%2 == 1)
      {
         this->mEven.bottomRows(1) *=0.5;
         this->mOdd.bottomRows(1).setZero();
      }
   }

   template <typename T> void InnerProduct<T>::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i)
   {
      rOut.row(i) = Precision::cast(ipolycol).transpose()*this->mIn;
   }

   template <typename T> void InnerProduct<T>::operator()(Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > rOut, const Eigen::Ref<const internal::Matrix>& ipolycol, const int i, const bool isEven)
   {
      if(isEven)
      {
         assert(ipolycol.rows() == this->mEven.rows());

         rOut.row(i) = Precision::cast(ipolycol).transpose()*this->mEven;
      } else
      {
         assert(ipolycol.rows() == this->mOdd.rows());

         rOut.row(i) = Precision::cast(ipolycol).transpose()*this->mOdd;
      }
   }
}
}
}
}

#endif // QUICC_POLYNOMIAL_LEGENDRE_EVALUATOR_INNERPRODUCT_HPP
