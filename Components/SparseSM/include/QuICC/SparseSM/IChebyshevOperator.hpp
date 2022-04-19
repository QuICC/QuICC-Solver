/** 
 * @file IChebyshevOperator.hpp
 * @brief Implementation of the generic interface to the Chebyshev sparse operator
 */

#ifndef QUICC_SPARSESM_ICHEBYSHEVOPERATOR_HPP
#define QUICC_SPARSESM_ICHEBYSHEVOPERATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/ISparseSMOperator.hpp"

namespace QuICC {

namespace SparseSM {

   /**
    * @brief Implementation of the generic interface to the Chebyshev sparse operator
    */ 
   class IChebyshevOperator: public ISparseSMOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IChebyshevOperator(const int rows, const int cols);

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevOperator();
         
      protected:
         /**
          * @brief Wrap around negative column indexes into the matrix
          */
         virtual void leftOutOfMatrix(TripletList_t& list, const int row, const int col, const Scalar_t value) const;

      private:
   };

}
}

#endif // QUICC_SPARSESM_ICHEBYSHEVOPERATOR_HPP
