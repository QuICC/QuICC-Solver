/**
 * @file GeostrophicSolidBody.cpp
 * @brief Source of the implementation of the solid body with unit angular momentum operator for a geostrophic basis
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "GeostrophicSolidBody.hpp"
#include "QuICC/Polynomial/Bessel/Generic.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "DenseSM/Bessel/details/GeostrophicTools.hpp"
#include "DenseSM/Bessel/GeostrophicAngularMomentum.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   GeostrophicSolidBody::GeostrophicSolidBody(const int nN, const Scalar_t sDNu)
      : IMatrixSMOperator(nN, 1), mNn(nN), mSDNu(sDNu)
   {
   }

   void GeostrophicSolidBody::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      assert(cols == 1);
      this->buildGenericOp(mat, rows);
   }

   void GeostrophicSolidBody::buildGenericOp(Internal::Matrix& mat, const int rows) const
   {
      mat = Internal::Matrix::Zero(rows, 1);

      Internal::Array isg, isw;
      details::GeostrophicTools::computeGridS(isg, isw, 2*rows, this->mSDNu);

      Internal::Matrix ipoly;
      ipoly.resize(isg.size(), rows);
      Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> jnl(this->mSDNu);
      jnl.compute<Internal::MHDFloat>(ipoly, rows, 1, isg, isw);

      mat.col(0) = (isg.transpose()*ipoly).transpose();

      // Normalize by angular momentum to get solid body with unit angular momentum
      GeostrophicAngularMomentum ag(rows, this->mSDNu);
      auto agOp = ag.mat();

      auto agMom = (mat.col(0).transpose()*agOp.col(0)).value();
      mat.col(0) /= agMom;
   }

} // Bessel
} // DenseSM
} // QuICC
