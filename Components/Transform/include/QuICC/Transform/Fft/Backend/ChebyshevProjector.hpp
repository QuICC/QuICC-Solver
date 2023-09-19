/**
 * @file ChebyshevProjector.hpp
 * @brief Interface for a generic API for Chebyshev projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"
#include "QuICC/Transform/Fft/Backend/StorageKind.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   /**
    * @brief Interface for a generic API Chebyshev projector
    *
    * Backward transform, spectral to physical space
    */
   class ChebyshevProjector
   {
      public:
         /// Typedef for the configuration class
         typedef Chebyshev::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Chebyshev::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         ChebyshevProjector();

         /**
          * @brief Destructor
          */
         ~ChebyshevProjector();

         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;

         /**
          * @brief Set Scaler array
          */
         void setScaler(const Array& scaler) const;

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
          * @brief Set output
          */
         void output(Matrix& rOut) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const;

         /**
          * @brief Set output scaled by array
          */
         void outputScale(Matrix& rOut) const;

         /**
          * @brief Set output scaled by array
          */
         void outputScale(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const;

         /**
          * @brief Apply FFT spectral to physical space
          */
         void applyFft(Matrix& phys, const Matrix& mods) const;

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
         void getSolution(Matrix& tmp, const int zeroRows = 0,
            const int extraRows = 0,
            const bool updateSolver = false) const;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_CHEBYSHEVPROJECTOR_HPP
