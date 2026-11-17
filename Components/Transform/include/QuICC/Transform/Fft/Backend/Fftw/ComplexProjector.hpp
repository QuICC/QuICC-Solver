/**
 * @file ComplexProjector.hpp
 * @brief Backend for a generic complex FFTW based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_COMPLEXPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_COMPLEXPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IComplexBackend.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Backend for a generic complex FFTW based projector
    */
   class ComplexProjector: public IComplexBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ComplexProjector() = default;

         /**
          * @brief Destructor
          */
         ~ComplexProjector() = default;

         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const override;

         /**
          * @brief Copy and pad FFT temporary input
          */
         void input(MatrixZ& tmp, const MatrixZ& in) const;

         /**
          * @brief Copy mean of field into backend
          */
         void inputMean(MatrixZ& tmp, const MatrixZ& in) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void inputDiff(MatrixZ& tmp, const MatrixZ& in, const int order, const MHDFloat scale) const;

         /**
          * @brief Scale with fast and slow index dependent functin
          */
         void inputDiff2D(MatrixZ& tmp, const MatrixZ& rData, const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks) const;

      protected:

      private:
         /**
          * @brief Force exact complex conjugate for mean
          */
         void forceConjugate(MatrixZ& rData) const;

         /**
          * @brief Apply padding
          */
         void applyPadding(MatrixZ& rData) const;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_COMPLEXPROJECTOR_HPP
