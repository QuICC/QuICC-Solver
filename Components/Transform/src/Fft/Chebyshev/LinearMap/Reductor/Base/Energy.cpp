/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Chebyshev energy reductor, with
 * linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Base/Energy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void Energy<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const
{
   this->mBackend.input(tmp, in);
}

void Energy<base_t>::applyPostOperator(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.output(rOut, tmp);
}

void Energy<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);

   // clean input from overflow errors.
   // Apparently not a problem in the Boussinesq case,
   // -> we do it for the anelastic case
   if(this->mBackend.getExtraSize() > 0)
   {
      tmp.bottomRows(tmp.rows()-this->mspSetup->specSize()).setZero();
   }   
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
