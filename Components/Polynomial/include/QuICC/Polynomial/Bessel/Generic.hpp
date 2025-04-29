/**
 * @file Generic.hpp
 * @brief Implementation of the spherical Bessel basis for genericl roots
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_GENERIC_HPP
#define QUICC_POLYNOMIAL_BESSEL_GENERIC_HPP

// System includes
//

// Project includes
//
#include "QuICC/Polynomial/Bessel/details/Operators.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

/**
 * @brief Implementation of the spherical Bessel basis for generic roots
 *
 * @tparam TOp Spherical bessel operator
 */
template <typename TOp> class Generic : private TOp
{
public:
   /**
    * @brief Default constructor
    *
    * @param dNu Bessel parameter: nu = l + dNu
    */
   Generic(const Internal::MHDFloat dNu);

   /**
    * @brief Compute spherical bessel basis for magnetic poloidal boundary
    * condition
    *
    * @param rOut    Output matrix
    * @param nPoly   Number of polynomials
    * @param l       Harmonic degree l
    * @param igrid   Physical grid points
    * @param scale   Scaling array, ie. weights
    */
   template <typename T>
   void compute(
      Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
      const int nPoly, const int l, const Internal::Array& igrid,
      const Internal::Array& scale);

protected:
   /**
    * @brief Bessel parameter nu = l + dNu
    */
   Internal::MHDFloat mDNu;
};

template <typename TOp>
Generic<TOp>::Generic(const Internal::MHDFloat dNu) : mDNu(dNu)
{}

template <typename TOp>
template <typename T>
inline void Generic<TOp>::compute(
   Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
   const int nPoly, const int l, const Internal::Array& igrid,
   const Internal::Array& scale)
{
   std::vector<Internal::MHDFloat> roots;
   details::getRoots(roots, l, nPoly + TOp::EXTRA_POLY, this->mDNu);

   if constexpr (TOp::EXTRA_L == 0)
   {
      TOp::compute(rOut, roots, l, igrid, scale, this->mDNu);
   }
   else
   {
      std::vector<Internal::MHDFloat> roots_extra;
      details::getRoots(roots_extra, l + TOp::EXTRA_L, nPoly + TOp::EXTRA_POLY,
         this->mDNu);

      TOp::compute(rOut, roots, roots_extra, l, igrid, scale, this->mDNu);
   }
}

} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_BESSEL_GENERIC_HPP
