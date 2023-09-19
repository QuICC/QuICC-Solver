/**
 * @file WorlandIntegrator.hpp
 * @brief Interface for a generic Worland cuFFT based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_WORLANDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_WORLANDINTEGRATOR_HPP

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
#include "QuICC/Transform/Fft/Backend/CuFft/IWorlandBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Worland cuFFT based integrator
    */
   class WorlandIntegrator: public IWorlandBackend
   {
      public:
         /**
          * @brief Constructor
          */
         WorlandIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~WorlandIntegrator();

         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity = false, const bool alwaysZeroNegative = false) const override;

         using IWorlandBackend::io;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(const bool isEven) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void input(const Matrix& in, const bool isEven) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void input(const MatrixZ& out, const bool isEven, const bool useReal) const;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut, const bool isEven) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const bool isEven, const bool useReal) const;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

         /**
          * @brief Convert FFT coefficients to Worland expansion
          */
         void forwardWorland(const bool isEven, const int id = 0) const;

         /**
          * @brief Lower alpha by 1
          */
         void lowerAlpha(const double alpha, const bool isEven, const int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower alpha by 1
          */
         void raiseAlpha(const double alpha, const bool isEven, const int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower beta by 1
          */
         void lowerBeta(const double alpha, const bool isEven, const int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Raise beta by 1
          */
         void raiseBeta(const double alpha, const bool isEven, const int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower beta by 1 and divide by r^2
          */
         void lowerR2Beta(const double alpha, const bool isEven, const int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Raise beta by 1 and multiply by r^2
          */
         void raiseR2Beta(const double alpha, const bool isEven, const int id = 0, const double norm = 1.0/std::sqrt(2.0), const bool scaleL0 = false) const;

         /**
          * @brief Apply sparse integration operator I2
          */
         void applyI2(const bool isEven, const int id = 0) const;

         /**
          * @brief Apply sparse integration operator I4
          */
         void applyI4(const bool isEven, const int id = 0) const;

      protected:
         /**
          * @brief Accurate expansion for given l
          */
         virtual int lSize(const int l) const override;

         /**
          * @brief Accurate expansion for given index i and l
          */
         virtual int lSize(const int i, const int l) const override;

      private:
         /**
          * @brief Apply even FFT
          */
         void applyEvenFft() const;

         /**
          * @brief Apply even FFT
          */
         void applyOddFft() const;

         /**
          * @brief Apply partial backward Worland transform
          */
         void partialForwardWorland(const int l, const int i0, const int start, const int cols, const GpuMatrix& J, double* out, const int oSize, const double normV, const double normM) const;

         /**
          * @brief FFT scaling factor
          */
         mutable double mFftScaling;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_WORLANDINTEGRATOR_HPP
