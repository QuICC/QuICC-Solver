/** 
 * @file MixedProjector.hpp
 * @brief Interface for a generic mixed FFTW based projector 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Backend/Fftw/IMixedBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic mixed FFTW based projector
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
          * @brief Initialise the FFTW transforms
          */
         virtual void init(const SetupType& setup) const override;

         /**
          * @brief Set input pointers for FFT
          */
         void input(const MatrixZ& in) const;

         /**
          * @brief Apply spectral derivative of given order
          */
         void inputDiff(const MatrixZ& rData, const int order, const MHDFloat scale) const;

         /**
          * @brief Set output data pointers for FFT
          */
         void output(MHDFloat* out) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MHDFloat* out, const MHDComplex* in) const;

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
         mutable MHDFloat* mpOut;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDPROJECTOR_HPP
