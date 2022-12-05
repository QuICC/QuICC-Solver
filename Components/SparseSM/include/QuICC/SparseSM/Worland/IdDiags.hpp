/**
 * @file IdDiags.hpp
 * @brief Interface to I2 diagonals for full sphere Worland (restricted) identity sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_IDDIAGS_HPP
#define QUICC_SPARSESM_WORLAND_IDDIAGS_HPP

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
    * @brief Implementation of the full sphere Worland (restricted) identity sparse operator
    */
   class IdDiags: public IDiags
   {
      public:
         /**
          * @brief Constructor
          */
         IdDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~IdDiags() = default;

         /**
          * @brief Main diagonal
          */
         virtual ACoeff_t d0(const ACoeff_t& n) const;

      protected:
         /**
          * @brief Restriction size
          */
         const int mQ;

      private:
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_IDDIAGS_HPP
