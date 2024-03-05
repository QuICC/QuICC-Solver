/**
 * @file r_1drSphJnlImplicit.hpp
 * @brief Implementation of the 1/r D r Bessel polynomial
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_R_1DRSPHJNLIMPLICIT_HPP
#define QUICC_POLYNOMIAL_BESSEL_R_1DRSPHJNLIMPLICIT_HPP

// System includes
//

// Project includes
//
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/Polynomial/Bessel/Tags.hpp"
#include "QuICC/Polynomial/Bessel/r_1drSphJnl.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

/// @brief Generic implementation
/// @tparam
template <class> class r_1drSphJnl;

/**
 * @brief Implementation of the Bessel polynomial with implicit grid division
 */
template <> class r_1drSphJnl<implicit_t>
{
public:
   /**
    * @brief Needed additional modes for computation
    */
   static const int EXTRA_POLY = 1;

   /**
    * @brief Additional roots with different harmonic degree
    */
   static const int EXTRA_L = -1;

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
      const std::vector<Internal::MHDFloat>& roots,
      const std::vector<Internal::MHDFloat>& roots_l1, const int l,
      const Internal::Array& igrid, const Internal::Array& scale,
      const Internal::MHDFloat dNu);
};

template <typename T>
inline void r_1drSphJnl<implicit_t>::compute(
   Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
   const std::vector<Internal::MHDFloat>& roots,
   const std::vector<Internal::MHDFloat>& roots_l1, const int l,
   const Internal::Array& igrid, const Internal::Array& scale,
   const Internal::MHDFloat dNu)
{
   Polynomial::Bessel::r_1drSphJnl<recurrence_t> r_1drSphJnl;

   if (l == 0)
   {
      // fallback on explicit implementation
      Polynomial::Bessel::r_1drSphJnl<explicit_t> r_1drSphJnlExplicit;
      r_1drSphJnlExplicit.compute<T>(rOut, roots, l, igrid, scale, dNu);
   }
   else
   {
      if (scale.size() == 0)
      {
         throw std::logic_error("this polynomial implementation needs to know "
                                "the weights for the internal projection");
      }

      const int nPoly = roots.size() - 1;

      if (igrid.size() < (2 * nPoly + l + 1) / 2)
      {
         throw std::logic_error(
            "Grid size does not allow exact integration for internal step");
      }

      assert(roots.size() == roots_l1.size());

      // Operates on polynomials with l = l-1
      int lm = std::abs(l - 1);
      int n_in = nPoly + 1;

      Internal::Matrix opA(igrid.size(), n_in);
      Polynomial::Bessel::SphJnl jnl;
      jnl.compute<Internal::MHDFloat>(opA, roots_l1, lm, igrid, scale, dNu);

      Internal::Matrix opB(igrid.size(), n_in);
      r_1drSphJnl.compute<Internal::MHDFloat>(opB, roots_l1, lm, igrid,
         Internal::Array(), dNu);

      Internal::Matrix opC(igrid.size(), nPoly);
      Polynomial::Bessel::SphJnl jnlB;
      std::vector<Internal::MHDFloat> ks(roots.begin(), std::prev(roots.end()));
      jnlB.compute<Internal::MHDFloat>(opC, ks, l, igrid, scale, dNu);

      rOut = ((opC.transpose() * opB * opA.transpose()).transpose()).cast<T>();
   }
}

} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_BESSEL_R_1DRSPHJNLIMPLICIT_HPP
