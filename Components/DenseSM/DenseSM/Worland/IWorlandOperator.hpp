/** 
 * @file IWorlandOperator.hpp
 * @brief Implementation of the generic interface to the full sphere Worland dense operator
 */

#ifndef QUICC_DENSESM_WORLAND_IWORLANDOPERATOR_HPP
#define QUICC_DENSESM_WORLAND_IWORLANDOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "DenseSM/IMatrixSMOperator.hpp"
#include "DenseSM/Worland/WorlandKind.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   /**
    * @brief Implementation of the generic interface to the full sphere Worland dense operator
    */ 
   class IWorlandOperator: public DenseSM::IMatrixSMOperator
   {
      public:
         /**
          * @brief Constructor with default Worland basis
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param alpha   Jacobi alpha parameter
          * @param dBeta   Jacobi dBeta parameter
          */
         IWorlandOperator(const int rows, const int cols);

         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param alpha   Jacobi alpha parameter
          * @param dBeta   Jacobi dBeta parameter
          */
         IWorlandOperator(const int rows, const int cols, const Scalar_t  alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~IWorlandOperator() = default;
         
      protected:
         /**
          * @brief Geostrophic alpha
          */
         const Scalar_t mcAlpha;

         /**
          * @brief Worland Jacobi beta = l + dBeta
          */
         const Scalar_t mcDBeta;

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

#endif // QUICC_DENSESM_WORLAND_IWORLANDOPERATOR_HPP
