/**
 * @file lowerSphJnl.hpp
 * @brief Implementation of the spherical Bessel basis
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_LOWERSPHJNL_HPP
#define QUICC_POLYNOMIAL_BESSEL_LOWERSPHJNL_HPP

// System includes
//

// Project includes
//
#include "QuICC/Polynomial/Bessel/details/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

/**
 * @brief Implementation of the spherical Bessel basis
 */
class lowerSphJnl
{
public:
   /**
    * @brief Needed additional modes for computation
    */
   static const int EXTRA_POLY = 0;

   /**
    * @brief Additional roots with different harmonic degree
    */
   static const int EXTRA_L = 0;

   /**
    * @brief Default constructor
    */
   lowerSphJnl() = default;

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
inline void lowerSphJnl::compute(
   Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
   const std::vector<Internal::MHDFloat>& roots, const int l,
   const Internal::Array& igrid, const Internal::Array& scale,
   const Internal::MHDFloat dNu)
{
   const int nPoly = roots.size();

   for (int j = 0; j < nPoly; j++)
   {
      auto k = roots.at(j);
      Internal::Array col(igrid.size());
      for (int i = 0; i < igrid.size(); i++)
      {
         col(i) = details::lowerSphJnl(k, l, igrid(i), dNu);
      }

      if (scale.size() > 0)
      {
         rOut.col(j).array() = (col.array() * scale.array()).cast<T>();
      }
      else
      {
         rOut.col(j).array() = col.cast<T>();
      }
   }
}

} // namespace Bessel
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_BESSEL_LOWERSPHJNL_HPP
