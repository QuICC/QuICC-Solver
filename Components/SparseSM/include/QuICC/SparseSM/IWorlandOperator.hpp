/** 
 * @file IWorlandOperator.hpp
 * @brief Implementation of the generic interface to the full sphere Worland sparse operator
 */

#ifndef QUICC_SPARSESM_IWORLANDOPERATOR_HPP
#define QUICC_SPARSESM_IWORLANDOPERATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/ISparseSMOperator.hpp"
#include "QuICC/SparseSM/Worland/WorlandKind.hpp"

namespace QuICC {

namespace SparseSM {

   /**
    * @brief Implementation of the generic interface to the full sphere Worland sparse operator
    */ 
   class IWorlandOperator: public ISparseSMOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandOperator(const int rows, const int cols, const Scalar_t  alpha, const Scalar_t dBeta);

         /**
          * @brief Destructor
          */
         virtual ~IWorlandOperator();
         
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

#endif // QUICC_SPARSESM_IWORLANDOPERATOR_HPP
