/** 
 * @file dPlm.hpp
 * @brief Implementation of the associated Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_DPLM_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_DPLM_HPP

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
    * @brief Implementation of the D associated Legendre polynomial
    */ 
   class dPlm: public ALegendreBase
   {
      public:
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);

         template <typename T, typename TEvaluator> void computedPl0(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEvaluator> void dPlm::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
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

      // Storage for P_l^{m-1} and P_l^{m+1}
      internal::Matrix iplm_1(gN, 2);
      internal::Matrix iplm1(gN, 2);

      internal::Matrix idiff(gN,1);
      if(m > 0)
      {
         // Initialize P_l^{m-1}
         ALegendreBase::Pmm(iplm_1.col(0), m-1, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            iplm_1.col(0).array() *= scale.segment(0,gN).array();
         }
         ALegendreBase::Pm1m(iplm_1.col(1), m-1, iplm_1.col(0), igrid, ALegendreBase::normPm1m());

         // Initialize \partial_theta P_l^m
         ALegendreBase::dPmm(idiff.col(0), m, iplm_1.col(1), ALegendreBase::normdPmm());
         evaluator(rOut, idiff.col(0), 0);

         if(nPoly > 1)
         {
            // Increment P_l^{m-1}
            ALegendreBase::Plm(iplm_1.col(0), m-1, m+1, iplm_1.col(1), iplm_1.col(0), igrid, ALegendreBase::normPlm());
            iplm_1.col(0).swap(iplm_1.col(1));

            // Increment P_l^{m+1}
            ALegendreBase::Pmm(iplm1.col(0), m+1, igrid, ALegendreBase::normPmm());
            if(scale.size() > 0)
            {
               iplm1.col(0).array() *= scale.segment(0,gN).array();
            }

            // Increment \partial_theta P_l^m
            ALegendreBase::dPlm(idiff.col(0), m, m+1, iplm_1.col(1), iplm1.col(0), ALegendreBase::normdPlm());
            evaluator(rOut, idiff.col(0), 1);
         }

         if(nPoly > 2)
         {
            // Increment P_l^{m-1}
            ALegendreBase::Plm(iplm_1.col(0), m-1, m+2, iplm_1.col(1), iplm_1.col(0), igrid, ALegendreBase::normPlm());
            iplm_1.col(0).swap(iplm_1.col(1));

            // Increment P_l^{m+1}
            ALegendreBase::Pm1m(iplm1.col(1), m+1, iplm1.col(0), igrid, ALegendreBase::normPm1m());

            // Increment \partial_theta P_l^m
            ALegendreBase::dPlm(idiff.col(0), m, m+2, iplm_1.col(1), iplm1.col(1), ALegendreBase::normdPlm());
            evaluator(rOut, idiff.col(0), 2);
         }

         for(int i = 3; i < nPoly; ++i)
         {
            int l = m + i;

            // Increment P_l^{m-1}
            ALegendreBase::Plm(iplm_1.col(0), m-1, l, iplm_1.col(1), iplm_1.col(0), igrid, ALegendreBase::normPlm());
            iplm_1.col(0).swap(iplm_1.col(1));

            // Increment P_l^{m+1}
            ALegendreBase::Plm(iplm1.col(0), m+1, l, iplm1.col(1), iplm1.col(0), igrid, ALegendreBase::normPlm());
            iplm1.col(0).swap(iplm1.col(1));

            // Increment \partial_theta P_l^m
            ALegendreBase::dPlm(idiff.col(0), m, l, iplm_1.col(1), iplm1.col(1), ALegendreBase::normdPlm());
            evaluator(rOut, idiff.col(0), i);
         }

      // m == 0 is a special case
      } else
      {
         this->computedPl0(rOut, nPoly, igrid, scale, evaluator);
      }
   }

   template <typename T, typename TEvaluator> void dPlm::computedPl0(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Storage for P_l^{m-1} and P_l^{m+1}
      internal::Matrix iplm1(gN, 2);
      internal::Matrix idiff(gN,1);

      // Initialize \partial_theta P_l^m
      idiff.col(0).setZero();
      evaluator(rOut, idiff.col(0), 0);

      if(nPoly > 1)
      {
         // Increment P_l^{m+1}
         ALegendreBase::Pmm(iplm1.col(0), 1, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            iplm1.col(0).array() *= scale.segment(0,gN).array();
         }

         // Increment \partial_theta P_l^m
         ALegendreBase::dPl0(idiff.col(0), 1, iplm1.col(0), ALegendreBase::normdPl0());
         evaluator(rOut, idiff.col(0), 1);
      }

      if(nPoly > 2)
      {
         // Increment P_l^{m+1}
         ALegendreBase::Pm1m(iplm1.col(1), 1, iplm1.col(0), igrid, ALegendreBase::normPm1m());

         // Increment \partial_theta P_l^m
         ALegendreBase::dPl0(idiff.col(0), 2, iplm1.col(1), ALegendreBase::normdPl0());
         evaluator(rOut, idiff.col(0), 2);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         // Increment P_l^{m+1}
         ALegendreBase::Plm(iplm1.col(0), 1, i, iplm1.col(1), iplm1.col(0), igrid, ALegendreBase::normPlm());
         iplm1.col(0).swap(iplm1.col(1));

         // Increment \partial_theta P_l^m
         ALegendreBase::dPl0(idiff.col(0), i, iplm1.col(1), ALegendreBase::normdPl0());
         evaluator(rOut, idiff.col(0), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_DPLM_HPP
