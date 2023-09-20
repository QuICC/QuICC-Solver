/**
 * @file MixedIntegrator.hpp
 * @brief Interface for a generic mixed cuFFT based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_MIXEDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_MIXEDINTEGRATOR_HPP

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
    * @brief Interface for a generic mixed cuFFT based integrator
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
          * @brief Initialise the cuFFT transforms
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
         void outputDiff(MatrixZ& rOut, const int order, const double scale) const;

         /**
          * @brief Apply derivative to ouput of given order, modified values can be specified for mixed operator
          */
         void outputDiff(MatrixZ& rOut, const int order, const double scale, const std::map<int,MHDComplex>& mod) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MatrixZ& out, const Matrix& in) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MHDComplex* out, const double* in) const;

      protected:

      private:
         /**
          * @brief Input data pointer
          */
         mutable const double* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable MHDComplex* mpOut;

         /**
          * @brief FFT scaling factor
          */
         mutable double mFftScaling;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_MIXEDINTEGRATOR_HPP
