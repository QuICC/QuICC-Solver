/** 
 * @file ChebyshevEnergy.hpp
 * @brief Interface for a generic Chebyshev cuFFT based energy reductor 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVENERGY_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVENERGY_HPP

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
#include "QuICC/Transform/Fft/Backend/CuFft/IChebyshevBackend.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Chebyshev cuFFT based energy reductor
    */ 
   class ChebyshevEnergy: public IChebyshevBackend
   {
      public:
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
         virtual void init(const SetupType& setup) const override;

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
         using IChebyshevBackend::io;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

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
          * @brief Apply padding
          */
         void applyPadding(Matrix& rData, const int extraRows = 0) const;

         /**
          * @brief Compute energy weights
          */
         void computeEWeights(const int size, const double lower, const double upper) const;

         /**
          * @brief Temporary data for mid operations
          */
         mutable Matrix  mTmpMid;

         /**
          * @brief Actual fft FWD size
          */
         mutable int mNrgFwdSize;

         /**
          * @brief Actual fft BWD size
          */
         mutable int mNrgBwdSize;

         /**
          * @brief Padding size
          */
         mutable int mPadSize;

         /**
          * @brief Solver for differential operators
          */
         mutable std::shared_ptr<Fftw::DifferentialSolver> mspSolver;

         /**
          * @brief FFT scaling factor
          */
         mutable double mFftScaling;

         /**
          * @brief Energy weights
          */
         mutable Array mEWeights;

         /**
          * @brief Spectral operator
          */
         mutable SparseMatrix mSpecOp;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_CHEBYSHEVENERGY_HPP
