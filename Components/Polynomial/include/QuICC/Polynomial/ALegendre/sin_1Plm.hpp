/**
 * @file sin_1Plm.hpp
 * @brief Implementation of the associated Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_SIN_1PLM_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_SIN_1PLM_HPP

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
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/ALegendre/ALegendreBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace ALegendre {

   /**
    * @brief Implementation of the associated Legendre polynomial
    */
   class sin_1Plm: public ALegendreBase
   {
      public:
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const Internal::Array& ifullgrid, const Internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void sin_1Plm::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const Internal::Array& ifullgrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      // Extract required part of grid
      int gN = (ifullgrid.rows()/2 + ifullgrid.rows()%2);
      Internal::Array igrid = ifullgrid.segment(0, gN);

      if (m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Polynomials is set to zero for m=0 as it only appears combined with \partial_\phi
      if(m == 0)
      {
         Internal::Matrix ipoly(gN, 1);
         ipoly.col(0).setZero();
         for(int i = 0; i < nPoly; ++i)
         {
            evaluator(rOut, ipoly.col(0), i, true);
         }

      } else
      {
         // Storage for P_{l+1}^{m+1} and P_{l+1}^{m-1}
         Internal::Matrix ipoly(gN,1);
         Internal::Matrix ipl1m1(gN,2);
         Internal::Matrix ipl1m_1(gN,2);

         // Initialize P_{l+1}^{m+1}
         ALegendreBase::Pmm(ipl1m1.col(0), m+1, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            ipl1m1.col(0).array() *= scale.segment(0,gN).array();
         }

         // Initialize P_{l+1}^{m-1}
         ALegendreBase::Pmm(ipl1m_1.col(0), m-1, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            ipl1m_1.col(0).array() *= scale.segment(0,gN).array();
         }
         ALegendreBase::Pm1m(ipl1m_1.col(1), m-1, ipl1m_1.col(0), igrid, ALegendreBase::normPm1m());
         ALegendreBase::Plm(ipl1m_1.col(0), m-1, m+1, ipl1m_1.col(1), ipl1m_1.col(0), igrid, ALegendreBase::normPlm());
         ipl1m_1.col(0).swap(ipl1m_1.col(1));

         // Initialize \frac{1}{sin(\theta)} P_{l}^{m}
         ALegendreBase::sin_1Plm(ipoly.col(0), m, m, ipl1m1.col(0), ipl1m_1.col(1), ALegendreBase::normsin_1Plm());
         evaluator(rOut, ipoly.col(0), 0, true);
         bool isEven = false;

         if(nPoly > 1)
         {
            // Increment P_{l+1}^{m+1}
            ALegendreBase::Pm1m(ipl1m1.col(1), m+1, ipl1m1.col(0), igrid, ALegendreBase::normPm1m());

            // Increment P_{l+1}^{m-1}
            ALegendreBase::Plm(ipl1m_1.col(0), m-1, m+2, ipl1m_1.col(1), ipl1m_1.col(0), igrid, ALegendreBase::normPlm());
            ipl1m_1.col(0).swap(ipl1m_1.col(1));

            // Increment \frac{1}{sin(\theta)} P_{l}^{m}
            ALegendreBase::sin_1Plm(ipoly.col(0), m, m+1, ipl1m1.col(1), ipl1m_1.col(1), ALegendreBase::normsin_1Plm());
            evaluator(rOut, ipoly.col(0), 1, isEven);
            isEven = !isEven;
         }

         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;

            // Increment P_{l+1}^{m+1}
            ALegendreBase::Plm(ipl1m1.col(0), m+1, l+1, ipl1m1.col(1), ipl1m1.col(0), igrid, ALegendreBase::normPlm());
            ipl1m1.col(0).swap(ipl1m1.col(1));

            // Increment P_{l+1}^{m-1}
            ALegendreBase::Plm(ipl1m_1.col(0), m-1, l+1, ipl1m_1.col(1), ipl1m_1.col(0), igrid, ALegendreBase::normPlm());
            ipl1m_1.col(0).swap(ipl1m_1.col(1));

            // Increment \frac{1}{sin(\theta)} P_{l}^{m}
            ALegendreBase::sin_1Plm(ipoly.col(0), m, l, ipl1m1.col(1), ipl1m_1.col(1), ALegendreBase::normsin_1Plm());
            evaluator(rOut, ipoly.col(0), i, isEven);
            isEven = !isEven;
         }
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_SIN_1PLM_HPP
