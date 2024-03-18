/**
 * @file CoriolisQp.cpp
 * @brief Source of the implementation of the full sphere Bessel Coriolis cross term acting on l-1
 */

// System includes
//
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "CoriolisQp.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"
#include "QuICC/Polynomial/Bessel/Generic.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Polynomial/Bessel/lowerSphJnl.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   CoriolisQp::CoriolisQp(const Internal::MHDFloat outDNu, const Internal::MHDFloat inDNu, const int rows, const int cols, const int l)
      : IEmbeddedSMOperator(rows, cols), mOutDNu(outDNu), mInDNu(inDNu), mL(l)
   {
   }

   void CoriolisQp::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      this->buildGenericOp(mat, rows, cols);
   }

   void CoriolisQp::buildGenericOp(Internal::Matrix& mat, const int rows, const int cols) const
   {
      const auto& l = this->mL;

      // Compute Legendre quadrature to integrate r polynomial
      int rp = 2*rows + l;
      int pts = 2*rp;
      Internal::Array igrid;
      Internal::Array iweights;
      Polynomial::Quadrature::WorlandSphEnergyRule wquad;
      wquad.computeQuadrature(igrid, iweights, pts);

      Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> jnl(this->mOutDNu);
      Polynomial::Bessel::Generic<Polynomial::Bessel::lowerSphJnl> qpJnl(this->mInDNu);

      Internal::Matrix tmpBwd(igrid.size(), rows);
      qpJnl.compute<Internal::MHDFloat>(tmpBwd, rows, l+1, igrid, Internal::Array());
      Internal::Matrix tmpFwd(igrid.size(), rows);
      jnl.compute<Internal::MHDFloat>(tmpFwd, rows, l, igrid, iweights.array());

      mat = tmpFwd.transpose()*tmpBwd;
   }

} // namespace Bessel
} // namespace DenseSM
} // namespace QuICC
