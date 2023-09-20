/**
 * @file Pn.hpp
 * @brief Implementation of the Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_LEGENDRE_PN_HPP
#define QUICC_POLYNOMIAL_LEGENDRE_PN_HPP

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
#include "Types/Precision.hpp"
#include "QuICC/Polynomial/Legendre/LegendreBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Legendre {

   /**
    * @brief Implementation of the Legendre polynomial
    */
   class Pn: public LegendreBase
   {
      public:
         /**
          * @brief Compute polynomial through recurrence relation
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::Array& fullGrid, const internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void Pn::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::Array& fullGrid, const internal::Array& scale, TEvaluator evaluator)
   {
      // Extract required part of grid
      int gN = (fullGrid.rows()/2 + fullGrid.rows()%2);
      internal::Array igrid = fullGrid.segment(0, gN);

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      internal::Matrix ipoly(gN, 2);

      LegendreBase::P0(ipoly.col(0), LegendreBase::normP0());
      if(scale.size() > 0)
      {
         ipoly.col(0).array() *= scale.segment(0,gN).array();
      }
      evaluator(rOut, ipoly.col(0), 0, true);
      bool isEven = false;

      if(nPoly > 1)
      {
         LegendreBase::P1(ipoly.col(1), igrid, LegendreBase::normP1());
         evaluator(rOut, ipoly.col(1), 1, isEven);
         isEven = !isEven;
      }

      for(int i = 2; i < nPoly; ++i)
      {
         LegendreBase::Pn(ipoly.col(0), i, ipoly.col(1), ipoly.col(0), igrid, LegendreBase::normPn());
         ipoly.col(0).swap(ipoly.col(1));
         evaluator(rOut, ipoly.col(1), i, isEven);
         isEven = !isEven;
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_LEGENDRE_PN_HPP
