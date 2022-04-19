/** 
 * @file MixedIntegrator.hpp
 * @brief Interface for a generic mixed FFTW based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDINTEGRATOR_HPP

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
    * @brief Interface for a generic mixed FFTW based integrator
    */ 
   class MixedIntegrator: public IMixedBackend
   {
      public:
         /**
          * @brief Constructor
          */
         MixedIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~MixedIntegrator();
         
         /**
          * @brief Initialise the FFTW transforms
          */
         virtual void init(const SetupType& setup) const override;
         
         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void output(MatrixZ& out) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale) const;

         /**
          * @brief Apply derivative to ouput of given order, modified values can be specified for mixed operator
          */
         void outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MatrixZ& out, const Matrix& in) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MHDComplex* out, const MHDFloat* in) const;

      protected:

      private:
         /**
          * @brief Input data pointer
          */
         mutable const MHDFloat* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable MHDComplex* mpOut;

         /**
          * @brief Bwd size
          */
         mutable int mBwdSize;

         /**
          * @brief Block size
          */
         mutable int mBlockSize;

         /**
          * @brief FFT scaling factor
          */
         mutable MHDFloat mFftScaling;

         /**
          * @brief Map for output data
          */
         mutable Eigen::Map<MatrixZ> mOutMap;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDINTEGRATOR_HPP
