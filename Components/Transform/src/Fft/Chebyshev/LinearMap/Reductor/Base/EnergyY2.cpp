/**
 * @file EnergyY2.cpp
 * @brief Source of the implementation of the Chebyshev energy Y^2 reductor,
 * with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Base/EnergyY2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void EnergyY2<base_t>::initOperator() const
{
   int extrasize = this->mBackend.getExtraSize(); // extrasize = 0 in the boussinesq case

   int inputCols = 2 * this->mspSetup->specSize() + extrasize;
   int size = inputCols + std::min(2, 2 * this->mspSetup->padSize());

   ::QuICC::SparseSM::Chebyshev::LinearMap::Y2 op(size, size,
      this->mspSetup->lower(), this->mspSetup->upper());

   this->mBackend.setSpectralOperator(
      op.mat().leftCols(inputCols));
}

void EnergyY2<base_t>::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in);
}

void EnergyY2<base_t>::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.outputSpectral(rOut, tmp);
}

void EnergyY2<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);

   // clean input from overflow errors.
   // Apparently not a problem in the Boussinesq case,
   // -> we do it for the anelastic case
   if(this->mBackend.getExtraSize() > 0)
   {
      tmp.bottomRows(tmp.rows() - this->mspSetup->specSize()).setZero();
   }

}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
