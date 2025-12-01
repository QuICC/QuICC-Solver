/**
 * @file d2Plm.hpp
 * @brief Implementation of the associated Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_ALEGENDRE_D2PLM_HPP
#define QUICC_POLYNOMIAL_ALEGENDRE_D2PLM_HPP

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
   class d2Plm: public ALegendreBase
   {
      public:
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

         template <typename T, typename TEvaluator> void computed2Pl0(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

         template <typename T, typename TEvaluator> void computed2Pl1(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEvaluator> void d2Plm::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int m, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      // Extract required part of grid
      int gN = igrid.rows();
      
      if (m < 0)
      {
         throw std::logic_error("Tried to compute associated Legendre polynomial derivative P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Storage for P_l^{m-2}, P_l^{m} and P_l^{m+2}
      Internal::Matrix iplm_2(gN, 2);
      Internal::Matrix iplm(gN, 2);
      Internal::Matrix iplm2(gN, 2);

      Internal::Matrix idiff(gN,1);

      if(m > 1)
      {
         // Initialize P_l^{m-2}
         ALegendreBase::Pmm(iplm_2.col(0), m-2, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            iplm_2.col(0).array() *= scale.segment(0,gN).array();
         }
         // Increment to P_m^{m-2}
         ALegendreBase::Pm1m(iplm_2.col(1), m-2, iplm_2.col(0), igrid, ALegendreBase::normPm1m());
         ALegendreBase::Plm(iplm_2.col(0), m-2, m, iplm_2.col(1), iplm_2.col(0), igrid, ALegendreBase::normPlm());
         iplm_2.col(0).swap(iplm_2.col(1));

         // Initialize P_l^{m}
         ALegendreBase::Pmm(iplm.col(0), m, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            iplm.col(0).array() *= scale.segment(0,gN).array();
         }

         // Initialize \partial2_theta P_l^m
         ALegendreBase::d2Pmm0m1(idiff.col(0), m, m, iplm_2.col(1), iplm.col(0), ALegendreBase::normd2Plm());
         evaluator(rOut, idiff.col(0), 0);

         if(nPoly > 1)
         {
            // Increment P_l^{m-2}
            ALegendreBase::Plm(iplm_2.col(0), m-2, m+1, iplm_2.col(1), iplm_2.col(0), igrid, ALegendreBase::normPlm());
            iplm_2.col(0).swap(iplm_2.col(1));

            // Increment P_l^{m}
            ALegendreBase::Pm1m(iplm.col(1), m, iplm.col(0), igrid, ALegendreBase::normPm1m());

            // Increment \partial2_theta P_l^m
            ALegendreBase::d2Pmm0m1(idiff.col(0), m, m+1, iplm_2.col(1), iplm.col(1), ALegendreBase::normd2Plm());
            evaluator(rOut, idiff.col(0), 1);

         }

         if(nPoly > 2)
         {
            // Increment P_l^{m-2}
            ALegendreBase::Plm(iplm_2.col(0), m-2, m+2, iplm_2.col(1), iplm_2.col(0), igrid, ALegendreBase::normPlm());
            iplm_2.col(0).swap(iplm_2.col(1));

            // Increment P_l^{m}
            ALegendreBase::Plm(iplm.col(0), m, m+2, iplm.col(1), iplm.col(0), igrid, ALegendreBase::normPlm());
            iplm.col(0).swap(iplm.col(1));

            // Initialize P_l^{m+2}
            ALegendreBase::Pmm(iplm2.col(0), m+2, igrid, ALegendreBase::normPmm());
            if(scale.size() > 0)
            {
               iplm2.col(0).array() *= scale.segment(0,gN).array();
            }

            // Increment \partial2_theta P_l^m
            ALegendreBase::d2Plm(idiff.col(0), m, m+2, iplm_2.col(1), iplm.col(1), iplm2.col(0), ALegendreBase::normd2Plm());
            evaluator(rOut, idiff.col(0), 2);
         }

         if(nPoly > 3)
         {
            // Increment P_l^{m-2}
            ALegendreBase::Plm(iplm_2.col(0), m-2, m+3, iplm_2.col(1), iplm_2.col(0), igrid, ALegendreBase::normPlm());
            iplm_2.col(0).swap(iplm_2.col(1));

            // Increment P_l^{m}
            ALegendreBase::Plm(iplm.col(0), m, m+3, iplm.col(1), iplm.col(0), igrid, ALegendreBase::normPlm());
            iplm.col(0).swap(iplm.col(1));

            // Increment P_l^{m+2}
            ALegendreBase::Pm1m(iplm2.col(1), m+2, iplm2.col(0), igrid, ALegendreBase::normPm1m());

            // Increment \partial2_theta P_l^m
            ALegendreBase::d2Plm(idiff.col(0), m, m+3, iplm_2.col(1), iplm.col(1), iplm2.col(1), ALegendreBase::normd2Plm());
            evaluator(rOut, idiff.col(0), 3);
         }

         for(int i = 4; i < nPoly; ++i)
         {
            int l = m + i;

            // Increment P_l^{m-2}
            ALegendreBase::Plm(iplm_2.col(0), m-2, l, iplm_2.col(1), iplm_2.col(0), igrid, ALegendreBase::normPlm());
            iplm_2.col(0).swap(iplm_2.col(1));

            // Increment P_l^{m}
            ALegendreBase::Plm(iplm.col(0), m, l, iplm.col(1), iplm.col(0), igrid, ALegendreBase::normPlm());
            iplm.col(0).swap(iplm.col(1));

            // Increment P_l^{m+2}
            ALegendreBase::Plm(iplm2.col(0), m+2, l, iplm2.col(1), iplm2.col(0), igrid, ALegendreBase::normPlm());
            iplm2.col(0).swap(iplm2.col(1));

            // Increment \partial2_theta P_l^m
            ALegendreBase::d2Plm(idiff.col(0), m, l, iplm_2.col(1), iplm.col(1), iplm2.col(1), ALegendreBase::normd2Plm());
            evaluator(rOut, idiff.col(0), i);

         }

      } 
      // m == 0 is special case
      else if(m == 0)
      {
         this->computed2Pl0(rOut, nPoly, igrid, scale, evaluator);
      }
      else if(m == 1)
      {
         this->computed2Pl1(rOut, nPoly, igrid, scale, evaluator);
      }
   }

   template <typename T, typename TEvaluator> void d2Plm::computed2Pl0(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      // Extract required part of grid
      int gN = igrid.rows();

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Storage for P_l^{0}, P_l^{2}
      Internal::Matrix ipl0(gN, 2);
      Internal::Matrix ipl2(gN, 2);

      Internal::Matrix idiff(gN,1);

      // Initialize P_l^{0}
      ALegendreBase::Pmm(ipl0.col(0), 0, igrid, ALegendreBase::normPmm());
      if(scale.size() > 0)
      {
         ipl0.col(0).array() *= scale.segment(0,gN).array();
      }

      // Initialize \partial2_theta P_l^0
      idiff.col(0).setZero();
      evaluator(rOut, idiff.col(0), 0);      

      if(nPoly > 1)
      {
         // Increment P_l^{0}
         ALegendreBase::Pm1m(ipl0.col(1), 0, ipl0.col(0), igrid, ALegendreBase::normPm1m());

         // Increment \partial2_theta P_l^0
         ALegendreBase::d2P10(idiff.col(0), ipl0.col(1), ALegendreBase::normd2P10());
         evaluator(rOut, idiff.col(0), 1);
      }

      if(nPoly > 2)
      {
         // Increment P_l^{0}
         ALegendreBase::Plm(ipl0.col(0), 0, 2, ipl0.col(1), ipl0.col(0), igrid, ALegendreBase::normPlm());
         ipl0.col(0).swap(ipl0.col(1));

         // Initialize P_l^{2}
         ALegendreBase::Pmm(ipl2.col(0), 2, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            ipl2.col(0).array() *= scale.segment(0,gN).array();
         }

         // Increment \partial2_theta P_l^0
         ALegendreBase::d2Pl0(idiff.col(0), 2, ipl2.col(0), ipl0.col(1), ALegendreBase::normd2Pl0());
         evaluator(rOut, idiff.col(0), 2);
         
      }

      if(nPoly > 3)
      {
         // Increment P_l^{0}
         ALegendreBase::Plm(ipl0.col(0), 0, 3, ipl0.col(1), ipl0.col(0), igrid, ALegendreBase::normPlm());
         ipl0.col(0).swap(ipl0.col(1));

         // Increment P_l^{2}
         ALegendreBase::Pm1m(ipl2.col(1), 2, ipl2.col(0), igrid, ALegendreBase::normPm1m());

         // Increment \partial2_theta P_l^0
         ALegendreBase::d2Pl0(idiff.col(0), 3, ipl2.col(1), ipl0.col(1), ALegendreBase::normd2Pl0());
         evaluator(rOut, idiff.col(0), 3);
         
      }

      for(int i = 4; i < nPoly; ++i)
      {
         // Increment P_l^{0}
         ALegendreBase::Plm(ipl0.col(0), 0, i, ipl0.col(1), ipl0.col(0), igrid, ALegendreBase::normPlm());
         ipl0.col(0).swap(ipl0.col(1));

         // Increment P_l^{2}
         ALegendreBase::Plm(ipl2.col(0), 2, i, ipl2.col(1), ipl2.col(0), igrid, ALegendreBase::normPlm());
         ipl2.col(0).swap(ipl2.col(1));

         // Increment \partial2_theta P_l^0
         ALegendreBase::d2Pl0(idiff.col(0), i, ipl2.col(1), ipl0.col(1), ALegendreBase::normd2Pl0());
         evaluator(rOut, idiff.col(0), i);
      }
   }

   template <typename T, typename TEvaluator> void d2Plm::computed2Pl1(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      // Extract required part of grid
      int gN = igrid.rows();
      
      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      // Storage for P_l^{1}, P_l^{3}
      Internal::Matrix ipl1(gN, 2);
      Internal::Matrix ipl3(gN, 2);

      Internal::Matrix idiff(gN,1);

      // Initialize P_l^{1}
      ALegendreBase::Pmm(ipl1.col(0), 1, igrid, ALegendreBase::normPmm());
      if(scale.size() > 0)
      {
         ipl1.col(0).array() *= scale.segment(0,gN).array();
      }
      // Initialize P_l^{3} to a zero value (for l<3)
      ipl3.col(0).setZero();
      ipl3.col(1).setZero();

      // Initialize \partial2_theta P_l^1
      ALegendreBase::d2Pl1(idiff.col(0), 1, ipl3.col(0), ipl1.col(0), ALegendreBase::normd2Pl1());
      evaluator(rOut, idiff.col(0), 0);

      if(nPoly > 1)
      {
         // Increment P_l^{1}
         ALegendreBase::Pm1m(ipl1.col(1), 1, ipl1.col(0), igrid, ALegendreBase::normPm1m());

         // Increment \partial2_theta P_l^1
         ALegendreBase::d2Pl1(idiff.col(0), 2, ipl3.col(0), ipl1.col(1), ALegendreBase::normd2Pl1()); 
         evaluator(rOut, idiff.col(0), 1);
      }

      if(nPoly > 2)
      {
         // Increment P_l^{1}
         ALegendreBase::Plm(ipl1.col(0), 1, 3, ipl1.col(1), ipl1.col(0), igrid, ALegendreBase::normPlm());
         ipl1.col(0).swap(ipl1.col(1));

         // Initialize P_l^{3}
         ALegendreBase::Pmm(ipl3.col(0), 3, igrid, ALegendreBase::normPmm());
         if(scale.size() > 0)
         {
            ipl3.col(0).array() *= scale.segment(0,gN).array();
         }

         // Increment \partial2_theta P_l^1
         ALegendreBase::d2Pl1(idiff.col(0), 3, ipl3.col(0), ipl1.col(1), ALegendreBase::normd2Pl1());
         evaluator(rOut, idiff.col(0), 2);
         
      }

      if(nPoly > 3)
      {
         // Increment P_l^{1}
         ALegendreBase::Plm(ipl1.col(0), 1, 4, ipl1.col(1), ipl1.col(0), igrid, ALegendreBase::normPlm());
         ipl1.col(0).swap(ipl1.col(1));

         // Increment P_l^{3}
         ALegendreBase::Pm1m(ipl3.col(1), 3, ipl3.col(0), igrid, ALegendreBase::normPm1m());

         // Increment \partial2_theta P_l^1
         ALegendreBase::d2Pl1(idiff.col(0), 4, ipl3.col(1), ipl1.col(1), ALegendreBase::normd2Pl1());
         evaluator(rOut, idiff.col(0), 3);
         
      }

      for(int i = 4; i < nPoly; ++i)
      {
         int l = 1 + i;

         // Increment P_l^{1}
         ALegendreBase::Plm(ipl1.col(0), 1, l, ipl1.col(1), ipl1.col(0), igrid, ALegendreBase::normPlm());
         ipl1.col(0).swap(ipl1.col(1));

         // Increment P_l^{3}
         ALegendreBase::Plm(ipl3.col(0), 3, l, ipl3.col(1), ipl3.col(0), igrid, ALegendreBase::normPlm());
         ipl3.col(0).swap(ipl3.col(1));

         // Increment \partial2_theta P_l^1
         ALegendreBase::d2Pl1(idiff.col(0), l, ipl3.col(1), ipl1.col(1), ALegendreBase::normd2Pl1());
         evaluator(rOut, idiff.col(0), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_ALEGENDRE_D2PLM_HPP
