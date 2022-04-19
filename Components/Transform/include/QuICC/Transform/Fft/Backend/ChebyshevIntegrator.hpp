/** 
 * @file ChebyshevIntegrator.hpp
 * @brief Interface for a generic API for Chebyshev integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVINTEGRATOR_HPP

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
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   /**
    * @brief Interface for a generic API for Chebyshev integrator
    */ 
   class ChebyshevIntegrator
   {
      public:
         /// Typedef for the configuration class
         typedef Chebyshev::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Chebyshev::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         ChebyshevIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevIntegrator();
         
         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void io(MHDFloat* out, const MHDFloat* in) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void io(Matrix& out, const Matrix& in) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void input(const MatrixZ& out, const bool useReal) const;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const bool useReal) const;

         /**
          * @brief Set output multiplied aby spectral operator
          */
         void outputSpectral(Matrix& rOut) const;

         /**
          * @brief Set output multiplied aby spectral operator
          */
         void outputSpectral(MatrixZ& rOut, const bool useReal) const;

         /**
          * @brief Apply FFT
          */
         void applyFft() const;

         /**
          * @brief Set spectral operator
          */
         void setSpectralOperator(const SparseMatrix& mat) const;

         /**
          * @brief Set spectral operator
          */
         void setMeanOperator(const SparseMatrix& mat) const;
         
      protected:

      private:
         /**
          * @brief PIMPL type forward declaration
          */
         struct BackendImpl;

         /**
          * @brief PIMPL
          */
         std::shared_ptr<BackendImpl> mpImpl;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVINTEGRATOR_HPP
