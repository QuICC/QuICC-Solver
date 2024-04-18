/**
 * @file Tor2Geostrophic.hpp
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */

#ifndef QUICC_DENSESM_BESSEL_TOR2GEOSTROPHIC_HPP
#define QUICC_DENSESM_BESSEL_TOR2GEOSTROPHIC_HPP

// System includes
//

// Project includes
//
#include "DenseSM/IMatrixSMOperator.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

/**
 * @brief Implementation of the projection operator from the toroidal scalar to
 * the geostrophic basis
 */
class Tor2Geostrophic : public IMatrixSMOperator
{
public:
   /**
    * @brief Constructor
    *
    * @param nN      Number of radial modes
    * @param nL      Number of harmonic degrees
    * @param nCpu    Number of CPU in MPI version
    * @param sDNu   Bessel dNu of S basis
    * @param torDNu Bessel dNu of toroidal basis
    */
   Tor2Geostrophic(const int nN, const int nL, const int nCpu, const Scalar_t sDNu,
      const Scalar_t torDNu);

   /**
    * @brief Destructor
    */
   virtual ~Tor2Geostrophic() = default;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    *
    * @param mat operator
    * @param rows rows of matrix
    * @param cols cols of matrix
    */
   void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const final;

   /**
    * @brief Max radial truncation
    */
   const int mNn;

   /**
    * @brief Number of harmonic degrees
    */
   const int mNl;

   /**
    * @brief Number of CPU
    */
   const int mNcpu;

   /**
    * @brief Bessel parameter nu = l + dnu for S grid
    */
   Internal::MHDFloat mSDNu;

   /**
    * @brief Bessel parameter nu = l + dnu for toroidal basis
    */
   Internal::MHDFloat mTorDNu;

private:
};

} // namespace Bessel
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_BESSEL_TOR2GEOSTROPHIC_HPP
