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
#include "QuICC/Transform/Fft/Backend/StorageKind.hpp"
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
          * @brief Copy input and pad
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          */
         void input(Matrix& tmp, const Matrix& in) const;

         /**
          * @brief Copy input real or imaginary part and shift
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          * @param shift
          */
         void input(Matrix& tmp, const Matrix& in, const int shift) const;

         /**
          * @brief Copy input real or imaginary part and pad
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          * @param useReal flag to extract real or im part
          */
         void input(Matrix& tmp, const MatrixZ& in, const bool useReal) const;

         /**
          * @brief Copy input real or imaginary part and shift
          *
          * @param tmp temporary storage
          * @param in input spectral coefficients
          * @param shift
          * @param useReal flag to extract real or im part
          */
         void input(Matrix& tmp, const MatrixZ& in, const int shift, const bool useReal) const;

         /**
          * @brief Square intermediate result
          */
         void square(Matrix& tmp, const Matrix& in, const bool isFirst) const;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut, const Matrix& tmp) const;

         /**
          * @brief Set output mutliplied by scalar operator
          */
         void outputSpectral(Matrix& rOut, const Matrix& tmp) const;

         /**
          * @brief Apply FFT
          */
         void applyFft(Matrix& phys, const Matrix& mods) const;

         /**
          * @brief Apply forward FFT
          */
         void applyFwdFft(Matrix& mods, const Matrix& phys) const;

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
         void getSolution(Matrix& tmp, const int zeroRows = 0, const int extraRows = 0) const;

         /**
          * @brief Get the temporary storage
          */
         Matrix& getStorage(const StorageKind = StorageKind::in) const;
         
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
