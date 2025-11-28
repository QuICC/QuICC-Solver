/**
 * @file Utils.cpp
 * @brief Source of the tools for Jones-Worland polynomial implementation
 */

// System includes
//

// Project includes
//
#include "QuICC/Polynomial/Worland/Utils.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

namespace Utils {

void selectJacobi(const std::string& t, Internal::MHDFloat& a, Internal::MHDFloat& db, Internal::Array& igrid, Internal::Array& iweights)
{
   auto defineWorland = [](const auto& wt, Internal::MHDFloat& a, Internal::MHDFloat& db, Internal::Array& igrid, Internal::Array& iweights)
   {
      // Set alpha and dBeta
      a = wt.ALPHA;
      db = wt.DBETA;

      // Compute quadrature
      int nR = igrid.size();
      if(nR > 0)
      {
         typename std::remove_reference<decltype(wt)>::type::Rule wquad;
         wquad.computeQuadrature(igrid, iweights, nR);
      }
   };

   if(t == "Chebyshev")
   {
      Polynomial::Worland::worland_chebyshev_t wt;
      defineWorland(wt, a, db, igrid, iweights);
   }
   else if(t == "Legendre")
   {
      Polynomial::Worland::worland_legendre_t wt;
      defineWorland(wt, a, db, igrid, iweights);
   }
   else if(t == "SphEnergy")
   {
      Polynomial::Worland::worland_sphenergy_t wt;
      defineWorland(wt, a, db, igrid, iweights);
   }
   else if(t == "CylEnergy")
   {
      Polynomial::Worland::worland_cylenergy_t wt;
      defineWorland(wt, a, db, igrid, iweights);
   }
   else
   {
      throw std::logic_error("Unknown Worland type");
   }
}

} // namespace Utils
} // namespace Worland
} // namespace Polynomial
} // namespace QuICC
