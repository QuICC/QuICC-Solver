/**
 * @file ChebyshevProjector.hpp
 * @brief Interface for a generic Chebyshev FFTW based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IChebyshevBackend.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic Chebyshev FFTW based projector
    * Backward transform, spectral to physical space
    */
   class ChebyshevProjector: public IChebyshevBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ChebyshevProjector() = default;

         /**
          * @brief Destructor
          */
         ~ChebyshevProjector() = default;

         /**
          * @brief Initialise the FFTW transforms
          */
         void init(const SetupType& setup) const final;

         /**
          * @brief Set Scaler array
          */
         void setScaler(const Array& scaler) const;

         /**
          * @brief Set output
          */
         void output(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& tmp, const bool useReal) const;

         /**
          * @brief Set output scaled by array
          */
         void outputScale(Eigen::Ref<Matrix> rOut) const;

         /**
          * @brief Set output scaled by array
          */
         void outputScale(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp, const bool useReal) const;

         /**
          * @brief Add new linear solver
          *
          * @param extra Extra spectral modes
          */
         void addSolver(const int extra = 0) const;

         /**
          * @brief Get solver
          */
         DifferentialSolver& solver() const;

         /**
          * @brief Get solution from solver
          *
          */
         void getSolution(Matrix& tmp, const int zeroRows = 0,
            const int extraRows = 0,
            const bool updateSolver = false) const;

         /**
          * @brief Apply padding
          */
         void applyPadding(Eigen::Ref<Matrix> rData, const int extraRows = 0) const final;
      protected:

      private:

         /**
          * @brief Solver for differential operators
          */
         mutable std::shared_ptr<DifferentialSolver> mspSolver;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_CHEBYSHEVPROJECTOR_HPP
