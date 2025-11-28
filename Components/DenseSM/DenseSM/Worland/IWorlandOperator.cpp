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
#include "include/QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

IWorlandOperator::IWorlandOperator(const int rows, const int cols,
   const Scalar_t alpha, const Scalar_t dBeta) :
    DenseSM::IMatrixSMOperator(rows, cols), mcAlpha(alpha), mcDBeta(dBeta)
{}

void IWorlandOperator::computeQuadrature(Internal::Array& igrid,
   Internal::Array& iweights, const int size) const
{
   // Using Chebyshev type
   {
      typedef Polynomial::Worland::worland_chebyshev_t WType;
      WType wt;
      if (this->mcAlpha == wt.ALPHA && this->mcDBeta == wt.DBETA)
      {
         WType::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
      }
   }

   // Using Legendre type
   {
      typedef Polynomial::Worland::worland_legendre_t WType;
      WType wt;
      if (this->mcAlpha == wt.ALPHA && this->mcDBeta == wt.DBETA)
      {
         WType::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
      }
   }

   // Using spherical energy type
   {
      typedef Polynomial::Worland::worland_sphenergy_t WType;
      WType wt;
      if (this->mcAlpha == wt.ALPHA && this->mcDBeta == wt.DBETA)
      {
         WType::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
      }
   }

   // Using cylindrical energy type
   {
      typedef Polynomial::Worland::worland_cylenergy_t WType;
      WType wt;
      if (this->mcAlpha == wt.ALPHA && this->mcDBeta == wt.DBETA)
      {
         WType::Rule quad;
         quad.computeQuadrature(igrid, iweights, size);
      }
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
