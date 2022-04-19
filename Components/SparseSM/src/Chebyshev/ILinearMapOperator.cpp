/** 
 * @file ILinearMapOperator.cpp
 * @brief Source of the implementation of generic interface to a chebyshev sparse operator using linear map y = ax + b, x = [-1, 1] (natural Chebyshev grid)
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

   ILinearMapOperator::ILinearMapOperator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper)
      : IChebyshevOperator(rows, cols), mA(-4242.0), mB(-4242.0)
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
