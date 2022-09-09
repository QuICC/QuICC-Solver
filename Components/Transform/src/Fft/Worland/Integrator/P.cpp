/**
 * @file P.cpp
 * @brief Source of the implementation of the Worland P integrator
 */

// System includes
//
#include <cassert>

// Debug includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Integrator/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   P::P()
   {
      this->setProfileTag();
   }

   P::~P()
   {
   }

   void P::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.forwardWorland(isEven);
   }

   void P::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven);
      this->mBackend.io(isEven);
   }

   void P::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->computeWorlandExpansion(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void P::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal);
      this->mBackend.io(isEven);
   }

   void P::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->computeWorlandExpansion(isEven);

      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
