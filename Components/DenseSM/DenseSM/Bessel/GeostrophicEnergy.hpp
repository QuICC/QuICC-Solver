/**
 * @file GeostrophicEnergy.hpp
 * @brief Implementation of the Energy operator for the geostrophic basis
 */

#ifndef QUICC_DENSESM_BESSEL_GEOSTROPHICENERGY_HPP
#define QUICC_DENSESM_BESSEL_GEOSTROPHICENERGY_HPP

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
    * @brief Implementation of the Energy operator for the geostrophic basis
    */
   class GeostrophicEnergy: public IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN      Number of radial modes
          * @param nL      Number of harmonic degrees
          * @param sDNu    Bessel dNu of S basis
          */
         GeostrophicEnergy(const int nN, const int nL, const Scalar_t sDNu);

         /**
          * @brief Destructor
          */
         virtual ~GeostrophicEnergy() = default;

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
          * @brief Bessel parameter nu = l + dnu for S grid
          */
         Internal::MHDFloat mSDNu;

      private:
   };

} // Bessel
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_BESSEL_GEOSTROPHICENERGY_HPP
