/** 
 * @file I6CylLaplhDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I6CylLaplhDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/I6CylLaplhDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   I6CylLaplhDiags::I6CylLaplhDiags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::I6CylLaplhDiags(alpha, MHD_MP(-0.5), l, q)
   {
      if(q > 0)
      {
         throw std::logic_error("I6CylLaplh: truncation for q>0 is not implemented");
      }
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d_5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 256.0*(l1 + n - 5.0).pow(2)*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)/((l1 + 2.0*n - 10.0)*(l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0));

      return this->normalizeDiag(n, -5)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d_4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = -128.0*(l1 + n - 4.0)*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)*(8.0*l2 + 4.0*l1*n - 33.0*l1 - 4.0*n.pow(2) + 16.0*n + 9.0)/((l1 + 2.0*n - 9.0)*(l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0));

      return this->normalizeDiag(n, -4)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d_3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      ACoeff_t val;

      val = 64.0*(l1 + n - 3.0)*(l1 + n - 2.0)*(l1 + n - 1.0)*(24.0*l3 - 4.0*l2*(6.0*n + 16.0) - 60.0*l1*n.pow(2) + 230.0*l1*n + 54.0*l1 - 12.0*n.pow(3) + 104.0*n.pow(2) - 183.0*n - 119.0)/((l1 + 2.0*n - 8.0)*(l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0));

      return this->normalizeDiag(n, -3)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = -128.0*(l1 + n - 2.0)*(l1 + n - 1.0)*(8.0*l4 - 40.0*l3*n - 10.0*l3 - 48.0*l2*n.pow(2) + 196.0*l2*n + 84.0*l2 + 16.0*l1*n.pow(3) + 52.0*l1*n.pow(2) - 236.0*l1*n - 157.0*l1 + 16.0*n.pow(4) - 64.0*n.pow(3) - 40.0*n.pow(2) + 208.0*n + 105.0)/((l1 + 2.0*n - 7.0)*(l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0));

      return this->normalizeDiag(n, -2)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      auto l5 = l3*l2;
      ACoeff_t val;

      val = 32.0*(l1 + n - 1.0)*(8.0*l5 - 152.0*l4*n - 24.0*l4 + 32.0*l3*n.pow(2) + 568.0*l3*n + 388.0*l3 + 448.0*l2*n.pow(3) - 272.0*l2*n.pow(2) - 1704.0*l2*n - 636.0*l2 + 272.0*l1*n.pow(4) - 944.0*l1*n.pow(3) - 832.0*l1*n.pow(2) + 2404.0*l1*n + 1299.0*l1 + 16.0*n.pow(5) - 240.0*n.pow(4) + 360.0*n.pow(3) + 800.0*n.pow(2) - 891.0*n - 585.0)/((l1 + 2.0*n - 6.0)*(l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0));

      return this->normalizeDiag(n, -1)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      auto l5 = l3*l2;
      ACoeff_t val;

      val = 16.0*(96.0*l5*n + 40.0*l5 - 480.0*l4*n.pow(2) - 800.0*l4*n - 360.0*l4 - 960.0*l3*n.pow(3) + 400.0*l3*n.pow(2) + 3600.0*l3*n + 1300.0*l3 + 2400.0*l2*n.pow(3) + 1920.0*l2*n.pow(2) - 4600.0*l2*n - 2880.0*l2 + 576.0*l1*n.pow(5) + 1200.0*l1*n.pow(4) - 3360.0*l1*n.pow(3) - 4600.0*l1*n.pow(2) + 3108.0*l1*n + 2035.0*l1 + 192.0*n.pow(6) - 1680.0*n.pow(4) + 3108.0*n.pow(2) - 675.0)/((l1 + 2.0*n - 5.0)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0));

      return this->normalizeDiag(n, 0)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      auto l4 = l2*l2;
      ACoeff_t val;

      val = 8.0*(2.0*n + 1.0)*(2.0*l1 + 2.0*n + 1.0)*(120.0*l4*n + 160.0*l4 - 160.0*l3*n.pow(2) - 760.0*l3*n - 900.0*l3 - 480.0*l2*n.pow(3) - 1120.0*l2*n.pow(2) + 1040.0*l2*n + 2240.0*l2 - 192.0*l1*n.pow(4) + 16.0*l1*n.pow(3) + 1912.0*l1*n.pow(2) + 804.0*l1*n - 2190.0*l1 + 16.0*n.pow(5) + 240.0*n.pow(4) + 360.0*n.pow(3) - 800.0*n.pow(2) - 891.0*n + 585.0)/((l1 + n)*(l1 + 2.0*n - 4.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0));

      return this->normalizeDiag(n, 1)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d2(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      auto l3 = l2*l1;
      ACoeff_t val;

      val = 8.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(40.0*l3*n + 90.0*l3 - 100.0*l2*n - 280.0*l2 - 48.0*l1*n.pow(3) - 244.0*l1*n.pow(2) - 156.0*l1*n + 365.0*l1 - 16.0*n.pow(4) - 64.0*n.pow(3) + 40.0*n.pow(2) + 208.0*n - 105.0)/((l1 + n)*(l1 + n + 1.0)*(l1 + 2.0*n - 3.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0));

      return this->normalizeDiag(n, 2)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d3(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = l1*l1;
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(60.0*l2*n + 190.0*l2 + 24.0*l1*n.pow(2) + 22.0*l1*n - 237.0*l1 - 12.0*n.pow(3) - 104.0*n.pow(2) - 183.0*n + 119.0)/((l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + 2.0*n - 2.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0));

      return this->normalizeDiag(n, 3)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d4(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(2.0*l1 + 2.0*n + 7.0)*(12.0*l1*n + 49.0*l1 + 4.0*n.pow(2) + 16.0*n - 9.0)/(2.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + n + 3.0)*(l1 + 2.0*n - 1.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0));

      return this->normalizeDiag(n, 4)*val;
   }

   I6CylLaplhDiags::ACoeff_t I6CylLaplhDiags::d5(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (n + 5.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*l1 + 2.0*n + 1.0)*(2.0*l1 + 2.0*n + 3.0)*(2.0*l1 + 2.0*n + 5.0)*(2.0*l1 + 2.0*n + 7.0)*(2.0*l1 + 2.0*n + 9.0)/(4.0*(l1 + n)*(l1 + n + 1.0)*(l1 + n + 2.0)*(l1 + n + 3.0)*(l1 + n + 4.0)*(l1 + 2.0*n + 1.0)*(l1 + 2.0*n + 2.0)*(l1 + 2.0*n + 3.0)*(l1 + 2.0*n + 4.0)*(l1 + 2.0*n + 5.0)*(l1 + 2.0*n + 6.0)*(l1 + 2.0*n + 7.0)*(l1 + 2.0*n + 8.0)*(l1 + 2.0*n + 9.0)*(l1 + 2.0*n + 10.0));

      return this->normalizeDiag(n, 5)*val;
   }

}
}
}
}
