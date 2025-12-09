/**
 * @file dsin_1Plm.hpp
 * @brief Implementation of the associated Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_DSIN_1PLM_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_DSIN_1PLM_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/ALegendre/ALegendreBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace ALegendre {

   /**
    * @brief Implementation of the D2 associated Legendre polynomial
    */
   class dsin_1Plm: public ALegendreBase
   {
      public:
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

         template <typename T, typename TEvaluator> void computedsin_1Pl1(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEvaluator> void dsin_1Plm::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();
      
      if (m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial derivative P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Storage for P_{l+1}^{m-2}, P_{l+1}^{m} and P_{l+1}^{m+2}
      Internal::Matrix ipl1m_2(gN, 2);
      Internal::Matrix ipl1m(gN, 2);
      Internal::Matrix ipl1m2(gN, 2);

      Internal::Matrix ipoly(gN,1);

      if(m > 1)
      {  // l=m; l+1 =m+1
         // Initialize P_{l+1}^{m-2}
         ALegendreBase::Pmm(ipl1m_2.col(0), m-2, igrid, ALegendreBase::normPmm()); // P_{m-2}^{m-2}
         if(scale.size() > 0)
         {
            ipl1m_2.col(0).array() *= scale.segment(0,gN).array();
         }
         // Increment to P_{m+1}^{m-2}
         ALegendreBase::Pm1m(ipl1m_2.col(1), m-2, ipl1m_2.col(0), igrid, ALegendreBase::normPm1m()); // P_{m-1}^{m-2}
         ALegendreBase::Plm(ipl1m_2.col(0), m-2, m, ipl1m_2.col(1), ipl1m_2.col(0), igrid, ALegendreBase::normPlm());  // P_{m}^{m-2}
         ipl1m_2.col(0).swap(ipl1m_2.col(1));
         ALegendreBase::Plm(ipl1m_2.col(0), m-2, m+1, ipl1m_2.col(1), ipl1m_2.col(0), igrid, ALegendreBase::normPlm()); // P_{m+1}^{m-2}
         ipl1m_2.col(0).swap(ipl1m_2.col(1));

         // Initialize P_{l+1}^{m}
         ALegendreBase::Pmm(ipl1m.col(0), m, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            ipl1m.col(0).array() *= scale.segment(0,gN).array();
         }
         ALegendreBase::Pm1m(ipl1m.col(1), m, ipl1m.col(0), igrid, ALegendreBase::normPm1m());

         // Initialize \partial_theta P_l^m/sin(theta)
         ALegendreBase::dsin_1Pmm(ipoly.col(0), m, ipl1m_2.col(1), ipl1m.col(1), ALegendreBase::normdsin_1Plm());
         evaluator(rOut, ipoly.col(0), 0);

         if(nPoly > 1)
         {
            // Increment P_{l+1}^{m-2}
            ALegendreBase::Plm(ipl1m_2.col(0), m-2, m+2, ipl1m_2.col(1), ipl1m_2.col(0), igrid, ALegendreBase::normPlm());
            ipl1m_2.col(0).swap(ipl1m_2.col(1));

            // Increment P_{l+1}^{m}
            ALegendreBase::Plm(ipl1m.col(0), m, m+2, ipl1m.col(1), ipl1m.col(0), igrid, ALegendreBase::normPlm());
            ipl1m.col(0).swap(ipl1m.col(1));

            // Initialize P_{l+1}^{m+2}
            ALegendreBase::Pmm(ipl1m2.col(0), m+2, igrid, ALegendreBase::normPmm());
            if(scale.size() > 0)
            {
               ipl1m2.col(0).array() *= scale.segment(0,gN).array();
            }

            // Increment \partial_theta P_l^m/sin(theta)
            ALegendreBase::dsin_1Plm(ipoly.col(0), m, m+1, ipl1m_2.col(1), ipl1m.col(1), ipl1m2.col(0), ALegendreBase::normdsin_1Plm());
            evaluator(rOut, ipoly.col(0), 1);

         }

         if(nPoly > 2)
         {
            // Increment P_{l+1}^{m-2}
            ALegendreBase::Plm(ipl1m_2.col(0), m-2, m+3, ipl1m_2.col(1), ipl1m_2.col(0), igrid, ALegendreBase::normPlm());
            ipl1m_2.col(0).swap(ipl1m_2.col(1));

            // Increment P_{l+1}^{m}
            ALegendreBase::Plm(ipl1m.col(0), m, m+3, ipl1m.col(1), ipl1m.col(0), igrid, ALegendreBase::normPlm());
            ipl1m.col(0).swap(ipl1m.col(1));

            // Initialize P_{l+1}^{m+2}
            ALegendreBase::Pm1m(ipl1m2.col(1), m+2, ipl1m2.col(0), igrid, ALegendreBase::normPm1m());

            // Increment \partial_theta P_l^m/sin(theta)
            ALegendreBase::dsin_1Plm(ipoly.col(0), m, m+2, ipl1m_2.col(1), ipl1m.col(1), ipl1m2.col(1), ALegendreBase::normdsin_1Plm());
            evaluator(rOut, ipoly.col(0), 2);
         }

         for(int i = 3; i < nPoly; ++i)
         {
            int l = m + i;

            // Increment P_{l+1}^{m-2}
            ALegendreBase::Plm(ipl1m_2.col(0), m-2, l+1, ipl1m_2.col(1), ipl1m_2.col(0), igrid, ALegendreBase::normPlm());
            ipl1m_2.col(0).swap(ipl1m_2.col(1));

            // Increment P_{l+1}^{m}
            ALegendreBase::Plm(ipl1m.col(0), m, l+1, ipl1m.col(1), ipl1m.col(0), igrid, ALegendreBase::normPlm());
            ipl1m.col(0).swap(ipl1m.col(1));

            // Initialize P_{l+1}^{m+2}
            ALegendreBase::Plm(ipl1m2.col(0), m+2, l+1, ipl1m2.col(1), ipl1m2.col(0), igrid, ALegendreBase::normPlm());
            ipl1m2.col(0).swap(ipl1m2.col(1));

            // Increment \partial_theta P_l^m/sin(theta)
            ALegendreBase::dsin_1Plm(ipoly.col(0), m, l, ipl1m_2.col(1), ipl1m.col(1), ipl1m2.col(1), ALegendreBase::normdsin_1Plm());
            evaluator(rOut, ipoly.col(0), i);

         }

      } 
      // Polynomials is set to zero for m=0 as it only appears combined with \partial_\phi
      else if(m == 0)
      {
         Internal::Matrix ipoly(gN, 1);
         ipoly.col(0).setZero();
         for(int i = 0; i < nPoly; ++i)
         {
            evaluator(rOut, ipoly.col(0), i, true);
         }
      }
      // m == 1 is special case
      else if(m == 1)
      {
         this->computedsin_1Pl1(rOut, nPoly, igrid, scale, evaluator);
      }
   }

   template <typename T, typename TEvaluator> void dsin_1Plm::computedsin_1Pl1(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();
      
      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Storage for P_{l+1}^{m} and P_{l+1}^{m+2}
      Internal::Matrix ipl11(gN, 2);
      Internal::Matrix ipl13(gN, 2);

      Internal::Matrix ipoly(gN,1);

      // Initialize P_{l+1}^{m}
      ALegendreBase::Pmm(ipl11.col(0), 1, igrid, ALegendreBase::normPmm());
      if(scale.size() > 0)
      {
         ipl11.col(0).array() *= scale.segment(0,gN).array();
      }
      // Increment to P_{m+1}^{m-2}
      ALegendreBase::Pm1m(ipl11.col(1), 1, ipl11.col(0), igrid, ALegendreBase::normPm1m());

      // Initialize \partial_theta P_l^m/sin(theta)
      ALegendreBase::dsin_1P11(ipoly.col(0), ipl11.col(1), ALegendreBase::normdsin_1Pl1());
      evaluator(rOut, ipoly.col(0), 0);

      if(nPoly > 1)
      {
         // Increment P_{l+1}^{m-2}
         ALegendreBase::Plm(ipl11.col(0), 1, 3, ipl11.col(1), ipl11.col(0), igrid, ALegendreBase::normPlm());
         ipl11.col(0).swap(ipl11.col(1));

         // Initialize P_{l+1}^{m}
         ALegendreBase::Pmm(ipl13.col(0), 3, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            ipl13.col(0).array() *= scale.segment(0,gN).array();
         }

         // Increment \partial_theta P_l^m/sin(theta)
         ALegendreBase::dsin_1Pl1(ipoly.col(0), 2, ipl13.col(0), ipl11.col(1), ALegendreBase::normdsin_1Pl1());
         evaluator(rOut, ipoly.col(0), 1);

      }

      if(nPoly > 2)
      {
         // Increment P_{l+1}^{m-2}
         ALegendreBase::Plm(ipl11.col(0), 1, 4, ipl11.col(1), ipl11.col(0), igrid, ALegendreBase::normPlm());
         ipl11.col(0).swap(ipl11.col(1));

         // Initialize P_{l+1}^{m}
         ALegendreBase::Pm1m(ipl13.col(1), 3, ipl13.col(0), igrid, ALegendreBase::normPm1m());

         // Increment \partial_theta P_l^m/sin(theta)
         ALegendreBase::dsin_1Pl1(ipoly.col(0), 3, ipl13.col(1), ipl11.col(1), ALegendreBase::normdsin_1Pl1());
         evaluator(rOut, ipoly.col(0), 2);
      }

      for(int i = 3; i < nPoly; ++i)
      { 
         int l = 1 + i;

         // Increment P_{l+1}^{m-2}
         ALegendreBase::Plm(ipl11.col(0), 1, l+1, ipl11.col(1), ipl11.col(0), igrid, ALegendreBase::normPlm());
         ipl11.col(0).swap(ipl11.col(1));

         // Increment P_{l+1}^{m}
         ALegendreBase::Plm(ipl13.col(0), 3, l+1, ipl13.col(1), ipl13.col(0), igrid, ALegendreBase::normPlm());
         ipl13.col(0).swap(ipl13.col(1));

         // Increment \partial_theta P_l^m/sin(theta)
         ALegendreBase::dsin_1Pl1(ipoly.col(0), l, ipl13.col(1), ipl11.col(1), ALegendreBase::normdsin_1Pl1());
         evaluator(rOut, ipoly.col(0), i);

   }

}
}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_DSIN_1PLM_HPP
