/**
 * @file Value.cpp
 * @brief Source of the implementation of boundary value
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   Value::Value(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : ICondition(alpha, dBeta, l, 0)
   {
   }

   Value::ACoeff_t Value::compute(const int maxN, const int k, const bool normalized)
   {
      ACoeff_t val;
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            val = this->valueChebyshev(maxN, k);
            break;
         case WorlandKind::LEGENDRE:
            val = this->valueLegendre(maxN, k);
            break;
         case WorlandKind::CYLENERGY:
            val = this->valueCylEnergy(maxN, k);
            break;
         case WorlandKind::SPHENERGY:
            val = this->valueSphEnergy(maxN, k);
            break;
      }

      ACoeff_t n = (ACoeffI::LinSpaced(maxN+1, 0, maxN)).cast<Scalar_t>();

      if(normalized)
      {
         val *= this->invnorm(n);
      }

      return val;
   }

   Value::ACoeff_t Value::valueChebyshev(const int maxN, const int k)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);
      val(0) = MHD_MP(1.0);

      for(int i = 1; i <= maxN; i++)
      {
         auto d2i = static_cast<Scalar_t>(2*i);
         auto d2k = static_cast<Scalar_t>(2*k);
         val(i) = val(i-1)*(d2i - MHD_MP(1.0) + d2k)/d2i;
      }

      return val;
   }

   Value::ACoeff_t Value::valueLegendre(const int maxN, const int k)
   {
      ACoeff_t val = ACoeff_t::Zero(maxN+1);

      if(k == 0)
      {
         val.setOnes();
      }
      else
      {
         val(0) = MHD_MP(1.0);

         for(int i = 1; i <= maxN; i++)
         {
            auto num = static_cast<Scalar_t>(i + k);
            auto den = static_cast<Scalar_t>(i);
            val(i) = val(i-1)*num/den;
         }
      }

      return val;
   }

   Value::ACoeff_t Value::valueCylEnergy(const int maxN, const int k)
   {
      // value depends only on alpha and n
      return this->valueLegendre(maxN, k);
   }

   Value::ACoeff_t Value::valueSphEnergy(const int maxN, const int k)
   {
      // value depends only on alpha and n
      return this->valueLegendre(maxN, k);
   }

} // Boundary
} // Worland
} // Polynomial
} // QuICC
