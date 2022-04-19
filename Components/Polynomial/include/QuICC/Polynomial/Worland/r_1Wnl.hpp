/** 
 * @file r_1Wnl.hpp
 * @brief Implementation of the 1/r Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_R_1WNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_R_1WNL_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the Worland polynomial
    */ 
   class r_1Wnl: public WorlandBase
   {
      public:
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void r_1Wnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland 1/r operator with l < 0");
      }

      if(nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if(gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      internal::Matrix ipoly(gN,2);

      internal::MHDFloat a = this->alpha(l);
      internal::MHDFloat b = this->beta(l);

      this->computeW0l(ipoly.col(0), l-1, a, b, igrid, WorlandBase::normWP0ab());
      if(scale.size() > 0)
      {
         ipoly.col(0).array() *= scale.array();
      }
      evaluator(rOut, ipoly.col(0), 0);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), a, b, ipoly.col(0), ixgrid, WorlandBase::normWP1ab());
         evaluator(rOut, ipoly.col(1), 1);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(0), i, a, b, ipoly.col(1), ipoly.col(0), ixgrid, WorlandBase::normWPnab());
         ipoly.col(0).swap(ipoly.col(1));
         evaluator(rOut, ipoly.col(1), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_R_1WNL_HPP
