/**
 * @file MixedIntegrator.hpp
 * @brief Interface for a generic mixed FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDINTEGRATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
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
         MixedIntegrator() = default;

         /**
          * @brief Destructor
          */
         ~MixedIntegrator() = default;

         /**
          * @brief Initialise the FFTW transforms
          */
         void init(const SetupType& setup) const override;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void output(Eigen::Ref<MatrixZ> out) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale) const;

         /**
          * @brief Apply derivative to ouput of given order, modified values can be specified for mixed operator
          */
         void outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const;

      protected:

      private:
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

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDINTEGRATOR_HPP
