/**
 * @file D1.cpp
 * @brief Source of the implementation of the ALegendre D1 integrator
 */

// External includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/ALegendre/Integrator/Base/D1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Integrator {

   D1<base_t>::D1()
   {
      this->setProfileTag();
   }


   void D1<base_t>::initBackend() const
   {
      int lshift = -1; // operator shifts l by one
      int extraN = 0; // no extra modes are required
      //this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void D1<base_t>::computeALegendreExpansion(const bool isEven) const
   {
      //this->mBackend.forwardALegendre(isEven);
      //this->mBackend.raiseR2Beta(-0.5, isEven);
   }

   void D1<base_t>::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      //this->mBackend.input(in, isEven);
      //this->mBackend.io(isEven);
   }

   void D1<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->computeALegendreExpansion(isEven);
      //this->mBackend.output(rOut, isEven);
   }

   void D1<base_t>::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      //this->mBackend.input(in, isEven, useReal);
      //this->mBackend.io(isEven);
   }

   void D1<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->computeALegendreExpansion(isEven);
      //this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
