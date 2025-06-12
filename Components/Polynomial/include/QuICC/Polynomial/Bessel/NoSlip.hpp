/**
 * @file NoSlip.hpp
 * @brief Implementation of the spherical Bessel basis for no-slip boundary
 * conditions. It consists in r^l and a spherical bessel basis.
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_NOSLIP_HPP
#define QUICC_POLYNOMIAL_BESSEL_NOSLIP_HPP

// System includes
//

// Project includes
//
#include "Polynomial/SphericalBessel/Jnl.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

/**
 * @brief Implementation of the spherical Bessel basis for no-slip boundary
 * conditions. It consists in r^l and a spherical bessel  basis
 */
template <typename TOp> class NoSlip : private TOp
{
public:
   /**
    * @brief Default constructor
    */
   NoSlip() = default;

   /**
    * @brief Compute spherical bessel basis for magnetic toroidal boundary
    * condition
    *
    * @param rOut    Output matrix
    * @param nPoly   Number of polynomials
    * @param l       Harmonic degree l
    * @param igrid   Physical grid points
    */
   template <typename T>
   void compute(
      Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
      const int nPoly, const int l, const Internal::Array& igrid,
      const Internal::Array& scale);
};

template <typename TOp>
template <typename T>
inline void NoSlip<TOp>::compute(
   Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
   const int nPoly, const int l, const Internal::Array& igrid,
   const Internal::Array& scale)
{
   std::vector<Internal::MHDFloat> roots = {0};
   SphericalBessel::getRoots(roots, l, nPoly - 1 + TOp::EXTRA_POLY, SphericalBessel::NoSlip_dNu());

   if constexpr (TOp::EXTRA_L == 0)
   {
      TOp::compute(rOut, roots, l, igrid, scale, SphericalBessel::NoSlip_dNu());
   }
   else
   {
      std::vector<Internal::MHDFloat> roots_extra = {};
      SphericalBessel::getRoots(roots_extra, l + TOp::EXTRA_L, nPoly + TOp::EXTRA_POLY,
         SphericalBessel::Value_dNu());

      TOp::compute(rOut, roots, roots_extra, l, igrid, scale,
         SphericalBessel::NoSlip_dNu(), SphericalBessel::Value_dNu());
   }
}

} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_BESSEL_NOSLIP_HPP
