/** 
 * @file ILinearMapOperator.hpp
 * @brief Implementation of the generic interface to the Cheyshev sparse operator based on a linear map y = ax + b, x = [-1, 1] (natural chebyshev grid)
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_ILINEARMAPOPERATOR_HPP
#define QUICC_SPARSESM_CHEBYSHEV_ILINEARMAPOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IChebyshevOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

   /**
    * @brief Implementation of the generic interface to the Chebyshev sparse operator based on a linear map y = ax + b, x = [-1, 1] (natural chebyshev grid)
    */ 
   class ILinearMapOperator: public IChebyshevOperator
   {
      public:
         /**
          * @brief Constructor
          * 
          * @param rows    Number of rows
          * @param cols    Number of columns 
          * @param lower   Lower bound of y
          * @param upper   Lower bound of y
          */
         ILinearMapOperator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~ILinearMapOperator();
         
      protected:
         /**
          * @brief Get mapping a coefficient from y = ax + b
          */
         Scalar_t a() const;

         /**
          * @brief Get mapping b coefficient from y = ax + b
          */
         Scalar_t b() const;

      private:
         /**
          * @brief Compute mapping
          *
          * @param lower   Lower bound
          * @param upper   Upper bound
          */
         void setBounds(const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief a coefficienct of y = ax + b
          */
         Scalar_t mA;

         /**
          * @brief b coefficienct of y = ax + b
          */
         Scalar_t mB;
   };

}
}
}

#endif // QUICC_SPARSESM_CHEBYSHEV_ILINEARMAPOPERATOR_HPP
