/**
 * @file MixedProjector.cpp
 * @brief Source of the interface for a generic FFTW based mixed projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/MixedProjector.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   MixedProjector::MixedProjector()
   {
   }

   MixedProjector::~MixedProjector()
   {
   }

   void MixedProjector::init(const SetupType& setup) const
   {
      //Initialize parent
      IMixedBackend::init(setup);

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();
      int blockSize = setup.blockSize();

      // Initialise temporary storage
      this->mTmp.setZero(bwdSize, blockSize);

      this->mPadSize = setup.padSize();

      // Create the two plans
      const int  *fftSize = &fwdSize;

      // create temporary storage for plan computation
      Matrix    tmpReal(fwdSize, blockSize);
      MatrixZ   tmpCplx(bwdSize, blockSize);

      // Create the complex to real plan
      this->mPlan = fftw_plan_many_dft_c2r(1, fftSize, blockSize, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, tmpReal.data(), NULL, 1, fwdSize, Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void MixedProjector::input(const MatrixZ& in) const
   {
      this->mTmp.topRows(this->mSpecSize) = in.topRows(this->mSpecSize);

      this->applyPadding(this->mTmp);
   }

   void MixedProjector::inputDiff(const MatrixZ& in, const int order, const MHDFloat scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*Math::cI;
         this->mTmp.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mSpecSize);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2);
         this->mTmp.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mSpecSize);
      }

      this->applyPadding(this->mTmp);
   }

   void MixedProjector::output(MHDFloat* out) const
   {
      this->io(out, this->mTmp.data());
   }

   void MixedProjector::io(MHDFloat* out, const MHDComplex* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

   void MixedProjector::applyPadding(MatrixZ& rData) const
   {
      // Set the m=0 values to zero
      rData.row(0).imag().setConstant(0);

      // Set the padded values to zero
      rData.bottomRows(this->mPadSize).setZero();
   }

   void MixedProjector::applyFft() const
   {
      fftw_execute_dft_c2r(this->mPlan, reinterpret_cast<fftw_complex* >(const_cast<MHDComplex *>(this->mpIn)), this->mpOut);
   }

}
}
}
}
}
