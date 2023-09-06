/**
 * @file IdDiags.hpp
 * @brief Interface to I2 diagonals for full sphere Worland (restricted) identity sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_IDDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_IDDIAGS_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland (restricted) identity sparse operator
    */
   class IdDiags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         IdDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~IdDiags() = default;

         /**
          * @brief Main diagonal
          *
          * @param n Array of n indexes
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const;

      protected:

      private:
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_IDDIAGS_HPP
