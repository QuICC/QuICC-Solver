/** 
 * @file MixedIntegrator.hpp
 * @brief Interface for a generic API for mixed integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_MIXEDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_MIXEDINTEGRATOR_HPP

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
    * @brief Interface for a generic API for mixed integrator
    */ 
   class MixedIntegrator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         MixedIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~MixedIntegrator();
         
         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;
         
         /**
          * @brief Apply FFT
          */
         void applyFft() const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void output(MatrixZ& out) const;

         /**
          * @brief Apply derivative to ouput of given order
          */
         void outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale) const;

         /**
          * @brief Apply derivative to ouput of given order, modified values can be specified for mixed operators
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

#endif // QUICC_TRANSFORM_FFT_BACKEND_MIXEDINTEGRATOR_HPP
