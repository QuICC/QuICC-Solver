/**
 * @file Operators.hpp
 * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_DETAILS_OPERATORS_HPP
#define QUICC_POLYNOMIAL_BESSEL_DETAILS_OPERATORS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

namespace details {

   /**
    * @brief Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat SphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief r of Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat rSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief 1/r of Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat r_1SphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief D of Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat dSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief D r of Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat drSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief 1/r D r of Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat r_1drSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief Spherical laplacian of Spherical Bessel basis Jnl(k, l, r)
    *
    * @param k Basis specific constant k
    * @param l Harmonic degree
    * @param r Radius r
    */
   Internal::MHDFloat slaplSphJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief Compute Bessel roots for Toroidal magnetic scalar
    */
   void getTorRoots(std::vector<Internal::MHDFloat>& roots, const int l, const int nRoots);

   /**
    * @brief Compute Bessel roots for Poloidal magnetic scalar
    */
   void getPolRoots(std::vector<Internal::MHDFloat>& roots, const int l, const int nRoots);

}
}
}
}

#endif // QUICC_POLYNOMIAL_BESSEL_DETAILS_OPERATORS_HPP
