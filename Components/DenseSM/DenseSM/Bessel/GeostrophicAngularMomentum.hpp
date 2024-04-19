/**
 * @file GeostrophicAngularMomentum.hpp
 * @brief Implementation of the angular momentum operator for the geostrophic basis
 */

#ifndef QUICC_DENSESM_BESSEL_GEOSTROPHICANGULARMOMENTUM_HPP
#define QUICC_DENSESM_BESSEL_GEOSTROPHICANGULARMOMENTUM_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/IMatrixSMOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   /**
    * @brief Implementation of the angula momentum operator for the geostrophic basis
    */
   class GeostrophicAngularMomentum: public IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN      Number of radial modes
          * @param sDNu    Bessel dNu of S basis
          */
         GeostrophicAngularMomentum(const int nN, const Scalar_t sDNu);

         /**
          * @brief Destructor
          */
         virtual ~GeostrophicAngularMomentum() = default;

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
          * @brief Bessel parameter nu = l + dnu for S grid
          */
         Internal::MHDFloat mSDNu;

      private:
         /**
          * @brief Build independent operator
          */
         void buildGenericOp(Internal::Matrix& mat, const int rows) const;
   };

} // Bessel
} // DenseSM
} // QuICC

#endif // QUICC_DENSESM_BESSEL_GEOSTROPHICANGULARMOMENTUM_HPP
