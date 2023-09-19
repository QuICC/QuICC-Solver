/**
 * @file ChebyshevProjector.hpp
 * @brief Interface for a generic Chebyshev cuFFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVPROJECTOR_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/IChebyshevBackend.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Chebyshev cuFFT based projector
    */
   class ChebyshevProjector: public IChebyshevBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ChebyshevProjector();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevProjector();

         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup) const override;

         /**
          * @brief Set Scaler array
          */
         void setScaler(const Array& scaler) const;

         /**
          * @brief Set input
          */
         void input(const Matrix& in, const bool needPadding = false) const;

         /**
          * @brief Set input
          */
         void input(const MatrixZ& in, const bool useReal, const bool needPadding = false) const;

         /**
          * @brief Set input and output to internal temporary storage
          */
         void io() const;
         using IChebyshevBackend::io;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const bool useReal) const;

         /**
          * @brief Set output scaled by array
          */
         void outputScale(Matrix& rOut) const;

         /**
          * @brief Set output scaled by array
          */
         void outputScale(MatrixZ& rOut, const bool useReal) const;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

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
         void getSolution(const int zeroRows = 0, const int extraRows = 0, const bool updateSolver = false) const;

      protected:

      private:
         /**
          * @brief Apply padding
          */
         void applyPadding(Matrix& rData, const int extraRows = 0) const;

         /**
          * @brief Padding size
          */
         mutable int mPadSize;

         /**
          * @brief Solver for differential operators
          */
         mutable std::shared_ptr<Fftw::DifferentialSolver> mspSolver;

         /**
          * @brief Scaler array
          */
         mutable Array mScaler;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVPROJECTOR_HPP
