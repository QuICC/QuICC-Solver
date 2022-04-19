/** 
 * @file IFftOperator.hpp
 * @brief Interface for a generic FFT based operator 
 */

#ifndef QUICC_TRANSFORM_FFT_IFFTOPERATOR_HPP
#define QUICC_TRANSFORM_FFT_IFFTOPERATOR_HPP

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
#include "QuICC/Transform/ITransformOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

   /**
    * @brief Interface for a generic FFT based operator
    */ 
   class IFftOperator: public ITransformOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IFftOperator();

         /**
          * @brief Destructor
          */
         virtual ~IFftOperator();

         /**
          * @brief Compute transform C2C or R2R componentwise
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const = 0;

         /**
          * @brief Compute transform R2C
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const = 0;

         /**
          * @brief Compute transform C2R
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const = 0;

         /**
          * @brief Compute transform R2R
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const = 0;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const = 0;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const = 0;

         /**
          * @brief Cleanup
          */
         void cleanup();

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;
         
      protected:
         /**
          * @brief Initialize base
          */
         void initBase() const;

      private:
         /**
          * @brief Operator specific initialization with empty default implementation
          */
         virtual void initOperator() const;

         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const = 0;
   };

}
}
}

#endif // QUICC_TRANSFORM_FFT_IFFTOPERATOR_HPP
