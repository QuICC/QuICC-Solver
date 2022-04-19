/** 
 * @file Plm.hpp
 * @brief Implementation of the associated Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_PLM_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_PLM_HPP

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
#include "QuICC/Polynomial/ALegendre/ALegendreBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace ALegendre {

   /**
    * @brief Implementation of the associated Legendre polynomial
    */ 
   class Plm: public ALegendreBase
   {
      public:
         /**
          * @brief Compute polynomial through recurrence relation
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const internal::Array& fullGrid, const internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void Plm::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const internal::Array& fullGrid, const internal::Array& scale, TEvaluator evaluator)
   {
      // Extract required part of grid
      int gN = (fullGrid.rows()/2 + fullGrid.rows()%2);
      internal::Array igrid = fullGrid.segment(0, gN);

      if (m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      internal::Matrix ipoly(gN, 2);

      ALegendreBase::Pmm(ipoly.col(0), m, igrid, ALegendreBase::normPmm());
      if(scale.size() > 0)
      {
         ipoly.col(0).array() *= scale.segment(0,gN).array();
      }
      evaluator(rOut, ipoly.col(0), 0, true);
      bool isEven = false;

      if(nPoly > 1)
      {
         ALegendreBase::Pm1m(ipoly.col(1), m, ipoly.col(0), igrid, ALegendreBase::normPm1m());
         evaluator(rOut, ipoly.col(1), 1, isEven);
         isEven = !isEven;
      }

      for(int i = 2; i < nPoly; ++i)
      {
         int l = m + i;
         ALegendreBase::Plm(ipoly.col(0), m, l, ipoly.col(1), ipoly.col(0), igrid, ALegendreBase::normPlm());
         ipoly.col(0).swap(ipoly.col(1));
         evaluator(rOut, ipoly.col(1), i, isEven);
         isEven = !isEven;
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_PLM_HPP
