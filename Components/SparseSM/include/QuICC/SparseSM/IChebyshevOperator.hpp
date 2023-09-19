/**
 * @file IChebyshevOperator.hpp
 * @brief Implementation of the generic interface to the Chebyshev sparse operator
 */

#ifndef QUICC_SPARSESM_ICHEBYSHEVOPERATOR_HPP
#define QUICC_SPARSESM_ICHEBYSHEVOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
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
          *
          * @param rows Number of rows
          * @param cols Number of columns
          */
         IChebyshevOperator(const int rows, const int cols);

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevOperator() = default;

      protected:
         /**
          * @brief Wrap around negative column indexes into the matrix
          *
          * @param list    List of Triplets representating the matrix
          * @param row     Row index
          * @param col     Column index
          * @param value   New value
          */
         virtual void leftOutOfMatrix(TripletList_t& list, const int row, const int col, const Scalar_t value) const;

      private:
   };

}
}

#endif // QUICC_SPARSESM_ICHEBYSHEVOPERATOR_HPP
