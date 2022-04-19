/**
 * @file MixedProjector.hpp
 * @brief Interface for a generic API for mixed projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_MIXEDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_MIXEDPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   /**
    * @brief Interface for a generic API for mixed projector
    */
   class MixedProjector
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         MixedProjector();

         /**
          * @brief Destructor
          */
         virtual ~MixedProjector();

         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;

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
         void applyFft() const;

         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief PIMPL type forward declaration
          */
         struct BackendImpl;

         /**
          * @brief PIMPL
          */
         std::shared_ptr<BackendImpl> mpImpl;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_MIXEDPROJECTOR_HPP
