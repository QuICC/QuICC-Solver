/** 
 * @file InsulatingShell.cpp
 * @brief Source of the implementation of the boundary value stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/InsulatingShell.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   InsulatingShell::InsulatingShell(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ISphericalOperator(rows, cols, lower, upper, l)
   {
   }

   InsulatingShell::ACoeff_t InsulatingShell::d_2(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto& l1 = this->l();
      const auto l2 = l1*l1;

      const auto val_num = a2*(2.0*l1*(l1+1.0)-2.0*(n-3.0)*n*((n-3.0)*n+5.0)-13.0)+a1*b1*(2.0*l1+1.0)*(2.0*(n-3.0)*n+5.0)+2.0*b2*(n.pow(2)-3.0*n+2.0).pow(2);
      const auto val_den = a2*(2.0*l1*(l1+1.0)-2.0*(n-1.0)*n*((n-1.0)*n+1.0)-1.0)+a1*b1*(2.0*l1+1.0)*(2.0*(n-1.0)*n+1.0)+2.0*b2*(n-1.0).pow(2)*n.pow(2);
      ACoeff_t val = -val_num/val_den;

      const auto corr_num = a1*(a1*(2.0*l2+2.0*l1-1.0)+2.0*b1*l1+b1);
      const auto corr_den = 2.0*(a2*(2.0*l2+2.0*l1-13.0)+5.0*a1*(2.0*b1*l1+b1)+8.0*b2);
      val(0) = -corr_num/corr_den;

      return val;
   }

   InsulatingShell::ACoeff_t InsulatingShell::d_1(const ACoeff_t& n) const
   {
      const auto& a1 = this->a();
      const auto a2 = a1*a1;
      const auto& b1 = this->b();
      const auto b2 = b1*b1;
      const auto& l1 = this->l();
      const auto l2 = l1*l1;

      const auto val_num = 4.0*a1*n*((2.0*l1+1.0)*a1 - b1);
      const auto val_den = a2*(2.0*l1*(l1+1.0)-2.0*n*(n+1.0)*(n.pow(2)+n+1.0)-1.0)+a1*b1*(2.0*l1+1.0)*(2.0*n*(n+1.0)+1.0)+2.0*b2*n.pow(2)*(n+1.0).pow(2);
      ACoeff_t val = val_num/val_den;

      const auto corr_num = 2.0*a1*(2.0*a1*l1+a1-b1);
      const auto corr_den = a2*(2.0*l2+2.0*l1-13.0)+5.0*a1*(2.0*b1*l1+b1)+8.0*b2;
      val(0) = corr_num/corr_den;

      return val;
   }

   InsulatingShell::ACoeff_t InsulatingShell::d0(const ACoeff_t& n) const
   {
      return ACoeff_t::Ones(n.size());
   }

   void InsulatingShell::buildTriplets(TripletList_t& list) const
   {
      ACoeffI ni = ACoeffI::LinSpaced(this->rows(), 0, this->rows()-1);
      ACoeff_t n = ni.cast<Scalar_t>();
      ACoeffI ni1 = ni.bottomRows(this->rows()-1);
      ACoeff_t n1 = n.bottomRows(this->rows()-1);
      ACoeffI ni2 = ni.bottomRows(this->rows()-2);
      ACoeff_t n2 = n.bottomRows(this->rows()-2);

      if(n.size() > 0)
      {
         list.reserve(3*std::max(this->rows(),this->cols()));
         this->convertToTriplets(list, -2, ni2, this->d_2(n2));
         this->convertToTriplets(list, -1, ni1, this->d_1(n1));
         this->convertToTriplets(list, 0, ni, this->d0(n));
      }
   }

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
