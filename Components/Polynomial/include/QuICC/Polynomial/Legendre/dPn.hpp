/**
 * @file dPn.hpp
 * @brief Implementation of the D Legendre polynomial
 */

#ifndef QUICC_POLYNOMIAL_LEGENDRE_DPN_HPP
#define QUICC_POLYNOMIAL_LEGENDRE_DPN_HPP

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
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Legendre/LegendreBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Legendre {

   /**
    * @brief Implementation of the D Legendre polynomial
    */
   class dPn: public LegendreBase
   {
      public:
         /**
          * @brief Compute polynomial through recurrence relation
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator);
   };

   template <typename T, typename TEvaluator> void dPn::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
   {
      int gN = igrid.rows();

      if (nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      Internal::Matrix ipn(gN, 2);
      Internal::Matrix idpn(gN, 2);

      idpn.col(0).setZero();
      evaluator(rOut, idpn.col(0), 0);

      if(nPoly > 1)
      {
         LegendreBase::P0(ipn.col(0), LegendreBase::normP0());
         if(scale.size() > 0)
         {
            ipn.col(0).array() *= scale.segment(0,gN).array();
         }

         LegendreBase::dP1(idpn.col(1), ipn.col(0), LegendreBase::normdP1());
         evaluator(rOut, idpn.col(1), 1);
      }

      if(nPoly > 2)
      {
         LegendreBase::P1(ipn.col(1), igrid, LegendreBase::normP1());

         LegendreBase::dPn(idpn.col(0), 2, ipn.col(1), idpn.col(1), igrid, LegendreBase::normdPn());
         idpn.col(0).swap(idpn.col(1));
         evaluator(rOut, idpn.col(1), 2);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         LegendreBase::Pn(ipn.col(0), i-1, ipn.col(1), ipn.col(0), igrid, LegendreBase::normPn());
         ipn.col(0).swap(ipn.col(1));

         LegendreBase::dPn(idpn.col(0), i, ipn.col(1), idpn.col(1), igrid, LegendreBase::normdPn());
         idpn.col(0).swap(idpn.col(1));
         evaluator(rOut, idpn.col(1), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_LEGENDRE_DPN_HPP
