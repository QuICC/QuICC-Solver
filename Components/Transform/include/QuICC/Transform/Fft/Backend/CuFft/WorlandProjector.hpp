/**
 * @file WorlandProjector.hpp
 * @brief Interface for a generic Worland cuFFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_WORLANDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_WORLANDPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Backend/CuFft/IWorlandBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Worland cuFFT based projector
    */
   class WorlandProjector: public IWorlandBackend
   {
      public:
         /**
          * @brief Constructor
          */
         WorlandProjector();

         /**
          * @brief Destructor
          */
         virtual ~WorlandProjector();

         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup, const int lshift, const bool lshiftOnlyParity = false, const bool alwaysZeroNegative = false) const override;

         /**
          * @brief Set input and output
          */
         void io(const bool isEven) const;

         /**
          * @brief Set input
          */
         void input(const Matrix& in, const bool isEven, const bool needPadding = false) const;

         /**
          * @brief Set input
          */
         void input(const MatrixZ& in, const bool isEven, const bool useReal, const bool needPadding = false) const;

         /**
          * @brief Set input and output to internal temporary storage
          */
         void io() const;
         using IWorlandBackend::io;

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
          * @brief Convert Worland expansion to  FFT coefficients
          */
         void backwardWorland(const bool isEven, const unsigned int id = 0) const;

         /**
          * @brief Lower alpha by 1
          */
         void lowerAlpha(const double alpha, const bool isEven, const unsigned int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower beta by 1
          */
         void lowerBeta(const double alpha, const bool isEven, const unsigned int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower beta by 1 and divide by r^2
          */
         void lowerR2Beta(const double alpha, const bool isEven, const unsigned int id = 0, const double norm = 1.0/std::sqrt(2.0)) const;

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
          * @brief Apply padding
          */
         void applyPadding(GpuMatrix& rData, const int extraRows = 0) const;

         /**
          * @brief Apply partial backward Worland transform
          */
         void partialBackwardWorland(const int l, const int i0, const int start, const int cols, const GpuMatrix& J, double* in, const int iSize, const double normV, const double normM) const;

         /**
          * @brief Padding size
          */
         mutable int mPadSize;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_WORLANDPROJECTOR_HPP
