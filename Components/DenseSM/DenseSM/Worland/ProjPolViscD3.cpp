/**
 * @file ProjPolViscD3.cpp
 * @brief Implementation of the poloidal projection of the term -D3 u_r/rho
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/ProjPolViscD3.hpp"
#include "DenseSM/Worland/DivR2F.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjPolViscD3::ProjPolViscD3(const int nNr, const int nNc, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{
   // bilaplacian
   this->mpOpA = std::make_shared<DivR2F>(nNr, nNc, lOut, mOut, lF, mF,
         lIn, mIn, pF, alpha, dBeta);
}

void ProjPolViscD3::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mpF->ls().size() != 1)
   {
      throw std::logic_error(
         "Operators are not implemented for forcing with multiple l");
   }

   const int Lout2 = this->mLout*(this->mLout+1);

   mat = -Lout2 * this->mpOpA->mat();
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
