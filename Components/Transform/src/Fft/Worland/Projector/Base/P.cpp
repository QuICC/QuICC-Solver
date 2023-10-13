/**
 * @file P.cpp
 * @brief Source of the implementation of the Worland P projector
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/P.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   P<base_t>::P()
   {
      this->setProfileTag();
   }

   void P<base_t>::computeWorlandExpansion(const bool isEven) const
   {
      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);
      this->mBackend.backwardWorland(isEven);
   }

   void P<base_t>::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void P<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void P<base_t>::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void P<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
