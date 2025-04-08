/**
 * @file EnergyD1.cpp
 * @brief Source of the implementation of the Chebyshev energy D^1 reductor,
 * with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Base/EnergyD1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void EnergyD1<base_t>::initOperator() const
{
   ::QuICC::SparseSM::Chebyshev::LinearMap::I1 op(
      this->mspSetup->specSize() + 1, this->mspSetup->specSize() + 1,
      this->mspSetup->lower(), this->mspSetup->upper());

   this->mBackend.solver().setOperator(op.mat());
}

void EnergyD1<base_t>::initBackend() const
{
   // Call parent initializer
   ILinearMapEnergy::initBackend();

   // Initialize the solver
   this->mBackend.addSolver();
}

void EnergyD1<base_t>::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in, 1);
   this->mBackend.getSolution(tmp, 1, 1);
}

void EnergyD1<base_t>::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.output(rOut, tmp);
}

void EnergyD1<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, 1, useReal);
   this->mBackend.getSolution(tmp, 1, 1);
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
