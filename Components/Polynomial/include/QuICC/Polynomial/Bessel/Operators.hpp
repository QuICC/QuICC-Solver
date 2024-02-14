/**
 * @file Operators.hpp
 * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_OPERATORS_HPP
#define QUICC_POLYNOMIAL_BESSEL_OPERATORS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

   /**
    * @brief Spherical Bessel basis Jnl(k, l, r)
    */
   Internal::MHDFloat Jnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief r of Spherical Bessel basis Jnl(k, l, r)
    */
   Internal::MHDFloat rJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief 1/r of Spherical Bessel basis Jnl(k, l, r)
    */
   Internal::MHDFloat r_1Jnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief D of Spherical Bessel basis Jnl(k, l, r)
    */
   Internal::MHDFloat dJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief D r of Spherical Bessel basis Jnl(k, l, r)
    */
   Internal::MHDFloat drJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

   /**
    * @brief Spherical laplacian of Spherical Bessel basis Jnl(k, l, r)
    */
   Internal::MHDFloat slaplJnl(const Internal::MHDFloat k, const int l, const Internal::MHDFloat r);

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

#endif // QUICC_POLYNOMIAL_BESSEL_OPERATORS_HPP
