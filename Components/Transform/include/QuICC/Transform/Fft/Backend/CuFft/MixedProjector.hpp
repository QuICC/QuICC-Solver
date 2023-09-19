/**
 * @file MixedProjector.hpp
 * @brief Interface for a generic mixed cuFFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_MIXEDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_MIXEDPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Backend/CuFft/IMixedBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic mixed cuFFT based projector
    */
   class MixedProjector: public IMixedBackend
   {
      public:
         /**
          * @brief Constructor
          */
         MixedProjector();

         /**
          * @brief Destructor
          */
         virtual ~MixedProjector();

         /**
          * @brief Initialise the cuFFT transforms
          */
         virtual void init(const SetupType& setup) const override;

         /**
          * @brief Set input pointers for FFT
          */
         void input(const MatrixZ& in) const;

         /**
          * @brief Apply spectral derivative of given order
          */
         void inputDiff(const MatrixZ& rData, const int order, const double scale) const;

         /**
          * @brief Set output data pointers for FFT
          */
         void output(double* out) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(double* out, const MHDComplex* in) const;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

      protected:

      private:
         /**
          * @brief Apply padding
          */
         void applyPadding(MatrixZ& rData) const;

         /**
          * @brief Padding size
          */
         mutable int mPadSize;

         /**
          * @brief Input data pointer
          */
         mutable const MHDComplex* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable double* mpOut;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_MIXEDPROJECTOR_HPP
