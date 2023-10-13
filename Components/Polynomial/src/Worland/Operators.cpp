/**
 * @file Tools.cpp
 * @brief Source of the tools for Jones-Worland polynomial implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Worland/Operators.hpp"

// Project includes
//
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Operators {

   void integrateRpWnl(Internal::Matrix& iop, const int l, const int p, const int nN)
   {
      int specSize = nN + l + p/2 + 1;
      int gridSize = specSize + 2;
      int lOut = l + p%2;

      // Create grids for integration
      Internal::Array igrid, iweights;
      Quadrature::LegendreRule lquad;
      lquad.computeQuadrature(igrid, iweights, gridSize);
      igrid.array() = ((igrid.array() + MHD_MP(1.0))/MHD_MP(2.0));

      // Compute integral weights
      Internal::Array volWeights(specSize,1);
      namespace ev = Evaluator;
      Wnl wnl;
      wnl.compute<Internal::MHDFloat>(volWeights, specSize, lOut, igrid, iweights, ev::Reduce());
      volWeights.array() *= MHD_MP(0.5);

      // Build operator
      Quadrature::WorlandRule quad;
      quad.computeQuadrature(igrid, iweights, gridSize);
      Internal::Matrix proj(gridSize, nN);
      wnl.compute<Internal::MHDFloat>(proj, nN, l, igrid, Internal::Array(), ev::Set());

      Internal::Matrix tmp(igrid.size(), specSize);
      wnl.compute<Internal::MHDFloat>(tmp, specSize, lOut, igrid, iweights, ev::Set());
      tmp.array().colwise() *= igrid.array().pow(p);
      iop = (volWeights.transpose()*tmp.transpose()*proj).cast<MHDFloat>().transpose();
   }

}
}
}
}
