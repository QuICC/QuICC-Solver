/** 
 * @file IEmbeddedOperator.hpp
 * @brief Implementation of the generic interface to the full sphere Worland dense operator
 */

#ifndef QUICC_DENSESM_WORLAND_IEMBEDDEDOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_IEMBEDDEDOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/IEmbeddedSMOperator.hpp"
#include "DenseSM/Worland/WorlandKind.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the generic interface to the full sphere Worland dense operator
    */ 
   class IEmbeddedOperator: public DenseSM::IEmbeddedSMOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param alpha   Jacobi alpha parameter
          * @param veta    Jacobi alpha parameter
          */
         IEmbeddedOperator(const int rows, const int cols, const Scalar_t  alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~IEmbeddedOperator() = default;
         
      protected:
         /**
          * @brief Type of Worland implementation
          */
         WorlandKind type() const;

      private:
         /**
          * Type of Worland implementation
          */
         WorlandKind mType;
   };

}
}
}

#endif // QUICC_DENSESM_WORLAND_IEMBEDDEDOPERATOR_HPP
