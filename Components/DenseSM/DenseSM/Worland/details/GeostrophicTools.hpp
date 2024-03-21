/**
 * @file GeostrophicTools.hpp
 * @brief Implementation of the base for geostrophic basis operator
 */

#ifndef QUICC_DENSESM_WORLAND_DETAILS_GEOSTROPHICTOOLS_HPP
#define QUICC_DENSESM_WORLAND_DETAILS_GEOSTROPHICTOOLS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

namespace details {

class GeostrophicTools
{
public:
   /**
    * @brief Akj coefficient used for geostrophic to toroidal projection
    * (Jiawen's thesis: (C.3))
    */
   static Internal::MHDFloat Akj(const int k, const int j);

   /**
    * @brief Bjn coefficient used for geostrophic to toroidal projection
    * (Jiawen's thesis: (3.28) but using Cn in place of Cnab)
    */
   static Internal::MHDFloat Bjn(const int j, const int n);

   /**
    * @brief Bjnab coefficient used for geostrophic to toroidal projection
    * (Jiawen's thesis: (3.28))
    */
   static Internal::MHDFloat Bjnab(const int j, const int n,
      const Internal::MHDFloat a, const Internal::MHDFloat b);

   /**
    * @Brief Normalization coefficient for $\Lambda_n(s)$ (Jiawen's thesis
    * (3.20))
    */
   static Internal::MHDFloat Cnab(const int n, const Internal::MHDFloat a,
      const Internal::MHDFloat b);

   /**
    * @Brief Normalization coefficient for $\tilde{\Lambda}_n(s)$ (Jiawen's
    * thesis (3.23))
    */
   static Internal::MHDFloat Cn(const int n);

   /**
    * @brief Compute quadrature points and weights for z integral
    */
   static void computeQuadratureZ(Internal::Array& igridz,
      Internal::Array& iweightz, const int nz);

   /**
    * @brief Compute Z integral
    */
   static void integrateZ(int l, int n, Internal::Matrix& iintgz,
      Internal::Array& igrids, Internal::Array& igridz,
      Internal::Array& iweightz);
};

} // namespace details
} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_DETAILS_GEOSTROPHICTOOLS_HPP
