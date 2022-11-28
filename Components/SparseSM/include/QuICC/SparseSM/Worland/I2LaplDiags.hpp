/**
 * @file I2LaplDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland I2Lapl sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I2LAPLDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_I2LAPLDIAGS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I2Lapl sparse operator
    */
   class I2LaplDiags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          */
         I2LaplDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I2LaplDiags() = default;

         /**
          * @brief 1. subdiagonal
          */
         virtual ACoeff_t d_1(const ACoeff_t& n) const = 0;

         /**
          * @brief Main diagonal
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const = 0;

         /**
          * @brief 1. superdiagonal
          */
         virtual ACoeff_t d1(const ACoeff_t& n) const = 0;

      protected:

      private:
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I2LAPLDIAGS_HPP
