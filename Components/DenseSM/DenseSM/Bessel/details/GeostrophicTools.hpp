/**
 * @file GeostrophicTools.hpp
 * @brief Implementation of the base for geostrophic basis operator
 */

#ifndef QUICC_DENSESM_BESSEL_DETAILS_GEOSTROPHICTOOLS_HPP
#define QUICC_DENSESM_BESSEL_DETAILS_GEOSTROPHICTOOLS_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

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
    * @brief Compute cylindrical S grid
    */
   static void computeGridS(Internal::Array& igridS, Internal::Array& iweightS,
      const int nS, const Internal::MHDFloat dNu);

   /**
    * @brief Compute quadrature points and weights for z integral
    */
   static void computeQuadratureZ(Internal::Array& igridZ,
      Internal::Array& iweightZ, const int nz);

   /**
    * @brief Compute Z integral
    *
    * @param l Harmonic degree
    * @param n Radial truncation
    * @param iintgz output z integral on s grid
    * @param igridS S grid
    * @param igridZ Z grid
    * @param iweightZ Z weights
    * @param dNu Bessel dNU
    */
   static void integrateZ(int l, int n, Internal::Matrix& iintgz,
      const Internal::Array& igridS, const Internal::Array& igridZ,
      const Internal::Array& iweightZ, const Internal::MHDFloat dNu);

   /**
    * @brief Cylindrical truncation nUg
    *
    * @param nN            Number of radial modes
    * @param nL            Number of harmonic degrees
    */
   static int cylTruncNug(const int nN, const int nL);

   /**
    * @brief Cylindrical truncation nS
    *
    * @param nr            Number of radial modes
    * @param nL            Number of harmonic degrees
    */
   static int cylTruncNs(const int nr, const int nL);

   /**
    * @brief Cylindrical truncation nZ
    *
    * @param nr            Number of radial modes
    * @param nL            Number of harmonic degrees
    */
   static int cylTruncNz(const int nr, const int nL);

   /**
    * @brief Cylindrical truncation nR
    *
    * @param nR            Number of radial modes
    * @param nL            Number of harmonic degrees
    */
   static int cylTruncNr(const int nR, const int nL);

   /**
    * @brief Cancel angular momentum
    *
    * @param spec       Field spectrum
    * @param momWeights Angular momentum weights per mode
    * @param solidBody  Solid body decomposition with unit angular momentum
    */
   static void cancelAngularMomentum(Array& spec, const Array& momWeights, const Array& solidBody);
};

} // namespace details
} // namespace Bessel
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_BESSEL_DETAILS_GEOSTROPHICTOOLS_HPP
