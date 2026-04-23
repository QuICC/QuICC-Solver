/**
 * @file Spectrum.cpp
 * @brief Source of the implementation of the Worland R^2 power spectrum operator
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/Base/Spectrum.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/InnerProduct.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/OuterProduct.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   Spectrum<base_t>::Spectrum()
      : IWorlandPower(1)
   {
      this->setProfileTag();
   }

   void Spectrum<base_t>::makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const
   {
      /*
      // copied from PowerR2
      int l = this->mspSetup->slow(i);
      int nPoly = this->mspSetup->fastSize(i);

      // Build operator
      op.resize(igrid.size(), nPoly);
      namespace ev = Polynomial::Worland::Evaluator;
      Polynomial::Worland::Wnl bwnl;
      bwnl.compute<MHDFloat>(op, nPoly, l, igrid, Internal::Array(), ev::Set());

      Polynomial::Worland::Wnl fwnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,Polynomial::Worland::worland_sphenergy_t::DBETA);

      eop.resize(igrid.size(), nPoly);
      fwnl.compute<MHDFloat>(eop, nPoly, l, igrid, iweights, ev::Set());
      */
      throw std::logic_error("Spectrum::makeOperator not implemented");
   }

   void Spectrum<base_t>::applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      // copied from PowerR2
      //this->defaultApplyOperator(rOut, i, in);
      throw std::logic_error("Spectrum::applyOperator not implemented");
   }

}
}
}
}
}
