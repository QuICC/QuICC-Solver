/** 
 * @file IWorlandOperator.hpp
 * @brief Implementation of the generic interface to the full sphere Worland dense operator
 */

#ifndef QUICC_DENSESM_IWORLANDOPERATOR_HPP
#define QUICC_DENSESM_IWORLANDOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "DenseSM/IDenseSMOperator.hpp"
#include "DenseSM/Worland/WorlandKind.hpp"

namespace QuICC {

namespace DenseSM {

   /**
    * @brief Implementation of the generic interface to the full sphere Worland dense operator
    */ 
   class IWorlandOperator: public IDenseSMOperator
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
         IWorlandOperator(const int rows, const int cols, const Scalar_t  alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~IWorlandOperator() = default;
         
      protected:
         /**
          * @brief Type of Worland implementation
          */
         Worland::WorlandKind type() const;

      private:
         /**
          * Type of Worland implementation
          */
         Worland::WorlandKind mType;
   };

}
}

#endif // QUICC_DENSESM_IWORLANDOPERATOR_HPP
