/**
 * @file ProjTorViscNu.cpp
 * @brief Implementation of the poloidal projection of the term \nu lapl(u) /rho
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/ProjTorViscNu.hpp"
#include "DenseSM/Worland/FC.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjTorViscNu::ProjTorViscNu(const int nNr, const int nNc, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{
   // bilaplacian
   this->mpOpA = std::make_shared<FC>(nNr, nNc, lOut, mOut, lF, mF,
         lIn, mIn, pF, alpha, dBeta);
}

void ProjTorViscNu::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mpF->ls().size() != 1)
   {
      throw std::logic_error(
         "Operators are not implemented for forcing with multiple l");
   }

   mat = this->mpOpA->mat();
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
