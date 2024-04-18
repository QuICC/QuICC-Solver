/**
 * @file Geostrophic2Tor.hpp
 * @brief Implementation of the projection operator from the geostrophic basis to Toroidal basis
 */

#ifndef QUICC_DENSESM_BESSEL_GEOSTROPHIC2TOR_HPP
#define QUICC_DENSESM_BESSEL_GEOSTROPHIC2TOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "DenseSM/IMatrixSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   /**
    * @brief Implementation of the projection operator from the geostrophic basis to Toroidal basis
    */
   class Geostrophic2Tor: public IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN      Number of radial modes
          * @param nL      Number of harmonic degrees
          * @param nCpu    Number of CPU in MPI version
          * @param sDNu    Bessel dNu of S basis
          * @param torDNu  Bessel dNu of toroidal basis
          */
         Geostrophic2Tor(const int nN, const int nL, const int nCpu, const Scalar_t sDNu, const Scalar_t torDNu);

         /**
          * @brief Destructor
          */
         virtual ~Geostrophic2Tor() = default;

      protected:
         /**
          * @brief Implementation of build dense matrix operator
          *
          * @param mat operator
          * @param rows rows of matrix
          * @param cols cols of matrix
          */
         void buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const final;

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

} // Bessel
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_BESSEL_GEOSTROPHIC2TOR_HPP
