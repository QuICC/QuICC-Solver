/**
 * @file Tor.hpp
 * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_TOR_HPP
#define QUICC_POLYNOMIAL_BESSEL_TOR_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Bessel/details/Operators.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

   /**
    * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
    */
   template <typename TOp> class Tor: private TOp
   {
      public:
         /**
          * @brief Default constructor
          */
         Tor() = default;

         /**
          * @brief Compute spherical bessel basis for magnetic toroidal boundary condition
          */
         template <typename T> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale);
   };

   template <typename TOp> template <typename T> inline void Tor<TOp>::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale)
   {
      std::vector<Internal::MHDFloat> roots;
      details::getTorRoots(roots, l, nPoly + TOp::EXTRA_POLY);

      if constexpr(TOp::EXTRA_L == 0)
      {
         TOp::compute(rOut, roots, l, igrid, scale);
      }
      else
      {
         std::vector<Internal::MHDFloat> roots_extra;
         details::getTorRoots(roots_extra, l + TOp::EXTRA_L, nPoly + TOp::EXTRA_POLY);

         TOp::compute(rOut, roots, roots_extra, l, igrid, scale);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_BESSEL_TOR_HPP
