/** 
 * @file ILinearMapOperator.cpp
 * @brief Source of the implementation of generic interface to a chebyshev sparse operator using linear map y = ax + b, x = [-1, 1] (natural Chebyshev grid)
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <limits>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

   ILinearMapOperator::ILinearMapOperator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : IChebyshevOperator(rows, cols), mA(std::numeric_limits<Scalar_t>::max()), mB(std::numeric_limits<Scalar_t>::max())
   {
      this->setBounds(lower,upper);
   }

   ILinearMapOperator::~ILinearMapOperator()
   {
   }

   void ILinearMapOperator::setBounds(const Scalar_t lower, const Scalar_t upper)
   {
      this->mA = (upper - lower)/MHD_MP(2.0);
      this->mB = (upper + lower)/MHD_MP(2.0);
   }

   IChebyshevOperator::Scalar_t ILinearMapOperator::a() const
   {
      return this->mA;
   }

   IChebyshevOperator::Scalar_t ILinearMapOperator::b() const
   {
      return this->mB;
   }

}
}
}
