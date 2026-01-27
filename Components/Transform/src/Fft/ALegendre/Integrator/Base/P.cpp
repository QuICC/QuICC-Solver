/**
 * @file P.cpp
 * @brief Source of the implementation of the ALegendre P integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/ALegendre/Integrator/Base/P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Integrator {

   P<base_t>::P()
   {
      this->setProfileTag();
   }


   void P<base_t>::computeALegendreExpansion(const bool isEven) const
   {
      this->mBackend.forwardALegendre(isEven);
   }

   void P<base_t>::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven);
      this->mBackend.io(isEven);
   }

   void P<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->computeALegendreExpansion(isEven);
      this->mBackend.output(rOut, isEven);
   }

   void P<base_t>::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal);
      this->mBackend.io(isEven);
   }

   void P<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->computeALegendreExpansion(isEven);

      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
