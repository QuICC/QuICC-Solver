/**
 * @file IWorlandOperator.cpp
 * @brief Source of the implementation of generic interface to a full sphere
 * Worland dense operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/IWorlandOperator.hpp"
#include "DenseSM/Worland/Tools.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

IWorlandOperator::IWorlandOperator(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta) :
    DenseSM::IMatrixSMOperator(rows, cols), mcAlpha(alpha), mcDBeta(dBeta)
{
   this->mType = Worland::Tools::identifyBasis(this->mcAlpha, this->mcDBeta);
}

void IWorlandOperator::computeQuadrature(Internal::Array& igrid,
   Internal::Array& iweights, const int size) const
{
   switch(this->type())
   {
      // Using Chebyshev type
      case WorlandKind::CHEBYSHEV:
      {
         Polynomial::Worland::worland_chebyshev_t::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
         break;
      }
      // Using Legendre type
      case WorlandKind::LEGENDRE:
      {
         Polynomial::Worland::worland_legendre_t::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
         break;
      }
      // Using cylindrical energy type
      case WorlandKind::CYLENERGY:
      {
         Polynomial::Worland::worland_cylenergy_t::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
         break;
      }
      // Using spherical energy type
      case WorlandKind::SPHENERGY:
      {
         Polynomial::Worland::worland_sphenergy_t::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
         break;
      }
   }
}

Worland::WorlandKind IWorlandOperator::type() const
{
   return this->mType;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
