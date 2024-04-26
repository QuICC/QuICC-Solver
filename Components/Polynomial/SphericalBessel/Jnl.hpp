/**
 * @file Jnl.hpp
 * @brief Implementation of the spherical Bessel
 * boundary conditions
 */

#ifndef QUICC_POLYNOMIAL_SPHERICALBESSEL_JNL_HPP
#define QUICC_POLYNOMIAL_SPHERICALBESSEL_JNL_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace SphericalBessel {

/**
 * @brief Zero boundary value Bessel function nu = l + dNu
 */
Internal::MHDFloat Value_dNu();

/**
 * @brief Insulating boundary Bessel function nu = l + dNu
 */
Internal::MHDFloat Insulating_dNu();

/**
 * @brief Norm of spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param dNu  Bessel nu = l + dNu
 */
Internal::MHDFloat norm(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat dNu);

/**
 * @brief Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief r of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat rSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief 1/r of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat r_1SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief D of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat dSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief D^2 of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat d2SphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief D r of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat drSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief 1/r D r of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat r_1drSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief Spherical laplacian of Spherical Bessel basis Jnl(k, l, r)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat slaplSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief Raise l by 1 on Spherical Bessel basis Jnl(k, l-1, r)
 *        [-D + l/r] f_l(z)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat raiseSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief Lower l by 1 on Spherical Bessel basis Jnl(k, l+1, r)
 *        [D + (l+1)/r] f_l(z)
 *
 * @param k Basis specific constant k
 * @param l Harmonic degree
 * @param r Radius r
 */
Internal::MHDFloat lowerSphJnl(const Internal::MHDFloat k, const int l,
   const Internal::MHDFloat r, const Internal::MHDFloat dNu);

/**
 * @brief Compute Bessel roots
 */
void getRoots(std::vector<Internal::MHDFloat>& roots, const int l,
   const int nRoots, const Internal::MHDFloat dNu);

} // namespace SphericalBessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_SPHERICALBESSEL_JNL_HPP
