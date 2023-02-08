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
#include "Profiler/Interface.hpp"

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
      this->mPlan = fftw_plan_many_dft_c2r(1, fftSize, blockSize, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, tmpReal.data(), NULL, 1, fwdSize, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void MixedProjector::input(MatrixZ& tmp, const MatrixZ& in) const
   {
      tmp.topRows(this->mSpecSize) = in.topRows(this->mSpecSize);

      this->applyPadding(tmp);
   }

   void MixedProjector::inputDiff(MatrixZ& out, const MatrixZ& in, const int order, const MHDFloat scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*Math::cI;
         out.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mSpecSize);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2);
         out.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mSpecSize);
      }

      this->applyPadding(out);
   }

   void MixedProjector::applyPadding(MatrixZ& rData) const
   {
      // Set the m=0 values to zero
      rData.row(0).imag().setConstant(0);

      // Set the padded values to zero
      rData.bottomRows(this->mPadSize).setZero();
   }

   void MixedProjector::applyFft(Matrix& phys, const MatrixZ& mods) const
   {
      Profiler::RegionFixture<4> fix("MixedProjector::applyFft");
      fftw_execute_dft_c2r(this->mPlan, reinterpret_cast<fftw_complex* >(const_cast<MHDComplex *>(mods.data())), phys.data());
   }

}
}
}
}
}
