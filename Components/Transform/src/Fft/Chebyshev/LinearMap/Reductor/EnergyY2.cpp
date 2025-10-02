/**
 * @file EnergyY2.cpp
 * @brief Source of the implementation of the Chebyshev energy Y^2 reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyY2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   EnergyY2::EnergyY2()
   {
   }

   EnergyY2::~EnergyY2()
   {
   }

   void EnergyY2::initOperator() const
   {
      // ExtraSize is = 0 in the Boussinesq case
      // ExtraSize is Background_Density::nN() in the anelastic case
      int extraSize = this->mBackend.getExtraSize(); 
      int size = 2*this->mspSetup->specSize() + extraSize + std::min(2, 2*this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::Y2 op(size, size, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setSpectralOperator(op.mat().leftCols(2*this->mspSetup->specSize() + extraSize));
   }

   void EnergyY2::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
   }

   void EnergyY2::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
   {
      assert(rOut.cols() == 1);
      this->mBackend.outputSpectral(rOut, tmp);
   }

   void EnergyY2::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      int extraSize = this->mBackend.getExtraSize();

      this->mBackend.input(tmp, in, useReal);
      
      // clean input from overflow errors.
      // Apparently not a problem in the Boussinesq case,
      // -> we do it for the anelastic case
      if(this->mBackend.getExtraSize() > 0)
      {
         tmp.bottomRows(tmp.rows()-this->mspSetup->specSize()).setZero();
      }   

   }

}
}
}
}
}
}
