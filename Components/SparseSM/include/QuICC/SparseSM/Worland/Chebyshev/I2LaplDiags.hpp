/**
 * @file I2LaplDiags.hpp
 * @brief Interface to I2Lapl diagonals for full sphere Worland I2Lapl sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_CHEBYSHEV_I2LAPLDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_CHEBYSHEV_I2LAPLDIAGS_HPP

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
#include "QuICC/SparseSM/Worland/I2LaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   /**
    * @brief Implementation of the full sphere Worland I2Lapl sparse operator
    */
   class I2LaplDiags: public QuICC::SparseSM::Worland::I2LaplDiags
   {
      public:
         /**
          * @brief Constructor
          */
         I2LaplDiags(const Scalar_t alpha, const int l);

         /**
          * @brief Destructor
          */
         virtual ~I2LaplDiags() = default;

         /**
          * @brief 1. subdiagonal
          */
         ACoeff_t d_1(const ACoeff_t& n) const final;

         /**
          * @brief Main diagonal
          */
         ACoeff_t d0(const ACoeff_t& n) const final;

         /**
          * @brief 1. superdiagonal
          */
         ACoeff_t d1(const ACoeff_t& n) const final;

      protected:

      private:
   };

}
}
}
}

#endif // QUICC_SPARSESM_WORLAND_CHEBYSHEV_I2LAPLDIAGS_HPP
