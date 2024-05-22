/**
 * @file GeostrophicSolidBody.hpp
 * @brief Implementation of the solid body with unit angular momentum decomposition 
 */

#ifndef QUICC_DENSESM_BESSEL_GEOSTROPHICSOLIDBODY_HPP
#define QUICC_DENSESM_BESSEL_GEOSTROPHICSOLIDBODY_HPP

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
    * @brief Implementation of the solid body with unit angular momentum decomposition
    */
   class GeostrophicSolidBody: public IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nN      Number of radial modes
          * @param sDNu    Bessel dNu of S basis
          */
         GeostrophicSolidBody(const int nN, const Scalar_t sDNu);

         /**
          * @brief Destructor
          */
         virtual ~GeostrophicSolidBody() = default;

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

#endif // QUICC_DENSESM_BESSEL_GEOSTROPHICSOLIDBODY_HPP
