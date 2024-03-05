/**
 * @file r_1drJnlExplicit.hpp
 * @brief Implementation of the 1/r D r Jnl
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_R_1DRSPHJNLEXPLICIT_HPP
#define QUICC_POLYNOMIAL_BESSEL_R_1DRSPHJNLEXPLICIT_HPP

// System includes
//

// Project includes
//
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Polynomial/Bessel/Tags.hpp"
#include "QuICC/Polynomial/Bessel/dSphJnl.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

/// @brief Generic implementation
/// @tparam
template <class> class r_1drSphJnl;

/**
 * @brief Implementation of the 1/r d r Jnl polynomial with explicit grid
 * division
 */
template <> class r_1drSphJnl<explicit_t>
{
public:
   /**
    * @brief Needed additional modes for computation
    */
   static const int EXTRA_POLY = 1;

   /**
    * @brief Additional roots with different harmonic degree
    */
   static const int EXTRA_L = 0;

   /**
    * @brief Compute spherical bessel basis
    *
    * @param rOut    Output matrix
    * @param roots   Vector of Bessel roots k: j(l, k_n r)
    * @param l       Harmonic degree l
    * @param igrid   Physical grid points
    * @param scale   Scaling array, ie. weights
    * @param dNu     nu shift used in bessel function J(l + dnu, x)
    */
   template <typename T>
   void compute(
      Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
      const std::vector<Internal::MHDFloat>& roots, const int l,
      const Internal::Array& igrid, const Internal::Array& scale,
      const Internal::MHDFloat dNu);
};

template <typename T>
inline void r_1drSphJnl<explicit_t>::compute(
   Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
   const std::vector<Internal::MHDFloat>& roots, const int l,
   const Internal::Array& igrid, const Internal::Array& scale,
   const Internal::MHDFloat dNu)
{
   if (scale.size() == 0)
   {
      throw std::logic_error("this polynomial implementation needs to know the "
                             "weights for the internal projection");
   }

   const int nPoly = roots.size() - 1;

   if (igrid.size() < (2 * nPoly + l + 2) / 2)
   {
      throw std::logic_error(
         "Grid size does not allow exact integration for internal step");
   }

   Internal::Matrix tOp(igrid.size(), nPoly);

   // Extend intermediate truncation by one due to multiplication by r
   Internal::Matrix opA(igrid.size(), nPoly + 1);
   Polynomial::Bessel::SphJnl jnl;
   jnl.compute<Internal::MHDFloat>(opA, roots, l, igrid,
      scale.array() * igrid.array(), dNu);

   Internal::Matrix opB(igrid.size(), nPoly + 1);
   Polynomial::Bessel::dSphJnl dJnl;
   dJnl.compute<Internal::MHDFloat>(opB, roots, l, igrid, Internal::Array(),
      dNu);

   Internal::Matrix opC(igrid.size(), nPoly);
   std::vector<Internal::MHDFloat> ks(roots.begin(), std::prev(roots.end()));
   jnl.compute<Internal::MHDFloat>(opC, ks, l, igrid,
      scale.array() * igrid.array().pow(-1), dNu);

   tOp = (opC.transpose() * opB * opA.transpose()).transpose();

   rOut = tOp.cast<T>();
}

} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_BESSEL_R_1DRSPHJNLEXPLICIT_HPP
