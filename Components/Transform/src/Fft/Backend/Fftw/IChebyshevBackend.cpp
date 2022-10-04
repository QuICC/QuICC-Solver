/**
 * @file IChebyshevBackend.cpp
 * @brief Source of the interface for a generic FFTW based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/IChebyshevBackend.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   IChebyshevBackend::IChebyshevBackend()
   {
   }

   IChebyshevBackend::~IChebyshevBackend()
   {
   }

   void IChebyshevBackend::init(const SetupType& setup) const
   {
      this->mSpecSize = setup.specSize();
      this->mBlockSize = setup.blockSize();
   }

   void IChebyshevBackend::input(Matrix& tmp, const Matrix& in) const
   {
      tmp.topRows(in.rows()) = in.topRows(in.rows());

      this->applyPadding(tmp);
   }

   void IChebyshevBackend::input(Matrix& tmp, const Matrix& in, const int shift) const
   {
      tmp.topRows(this->mSpecSize - shift) = in.block(shift, 0, this->mSpecSize - shift, in.cols());
   }

   void IChebyshevBackend::input(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const
   {
      if(useReal)
      {
         tmp.topRows(in.rows()) = in.real().topRows(in.rows());
      } else
      {
         tmp.topRows(in.rows()) = in.imag().topRows(in.rows());
      }

      this->applyPadding(tmp);
   }

   void IChebyshevBackend::input(Matrix& tmp, const MatrixZ& in,
      const int shift, const bool useReal) const
   {
      if(useReal)
      {
         tmp.topRows(this->mSpecSize - shift) = in.real().block(shift, 0, this->mSpecSize - shift, in.cols());
      } else
      {
         tmp.topRows(this->mSpecSize - shift) = in.imag().block(shift, 0, this->mSpecSize - shift, in.cols());
      }
      tmp.bottomRows(this->mSpecSize - shift).setZero();
   }

   void IChebyshevBackend::applyPadding(Matrix& rData, const int extraRows) const
   {
      if(this->mPadSize > 0)
      {
         throw std::logic_error("if needed it should be implemented in a derived class");
      }
   }

   Matrix& IChebyshevBackend::getStorage(const StorageKind kind) const
   {
      if (kind == StorageKind::in)
      {
         return this->mTmp;
      }
      else
      {
         return this->mTmpComp;
      }
   }

}
}
}
}
}
