/** 
 * @file IChebyshevEnergy.hpp
 * @brief Interface for a generic Chebyshev FFT based energy reductor 
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVENERGY_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVENERGY_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"
#include "QuICC/Transform/Fft/Backend/ChebyshevEnergy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   /**
    * @brief Interface for a generic Chebyshev FFT based energy reductor
    */ 
   class IChebyshevEnergy: public IChebyshevOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IChebyshevEnergy();

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevEnergy();

         /**
          * @brief Compute reduction of complex data
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute reduction of real data
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const override;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;
         
      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief FFT backend
          */
         Backend::ChebyshevEnergy mBackend;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& tmp, const Matrix& in) const = 0;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut, const Matrix& tmp) const = 0;

         /**
          * @brief Apply pre FFT operator for component wise operations
          *
          * @param in   Input values
          */
         virtual void applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const = 0;

         /**
          * @brief Compute transform R2C (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;

         /**
          * @brief Compute transform C2R (disabled)
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;
   };

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ICHEBYSHEVENERGY_HPP
