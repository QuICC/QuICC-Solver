/** 
 * @file ChebyshevEnergy.hpp
 * @brief Interface for a generic API for Chebyshev energy reductor 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVENERGY_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVENERGY_HPP

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
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   /**
    * @brief Interface for a generic API for Chebyshev energy reductor
    */ 
   class ChebyshevEnergy
   {
      public:
         /// Typedef for the configuration class
         typedef Chebyshev::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Chebyshev::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         ChebyshevEnergy();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevEnergy();
         
         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;

         /**
          * @brief set spectral operator
          */
         void setSpectralOperator(const SparseMatrix& mat) const;

         /**
          * @brief Set input
          */
         void input(const Matrix& in, const bool needPadding = false) const;

         /**
          * @brief Set input
          */
         void input(const MatrixZ& in, const bool useReal, const bool needPadding = false) const;

         /**
          * @brief Square intermediate result
          */
         void square(const bool isFirst) const;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut) const;

         /**
          * @brief Set output mutliplied by scalar operator
          */
         void outputSpectral(Matrix& rOut) const;

         /**
          * @brief Set input and output to internal temporary storage
          */
         void io() const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void io(MHDFloat* out, const MHDFloat* in) const;

         /**
          * @brief Apply FFT
          */
         void applyFft() const;

         /**
          * @brief Apply forward FFT
          */
         void applyFwdFft() const;

         /**
          * @brief Add new linear solver
          *
          * @param extra Extra spectral modes
          */
         void addSolver(const int extra = 0) const;

         /**
          * @brief Get solver
          */
         Fftw::DifferentialSolver& solver() const;

         /**
          * @brief Get solution from solver
          *
          */
         void getSolution(const int zeroRows = 0, const int extraRows = 0) const;
         
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

#endif // QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVENERGY_HPP
