/**
 * @file MixedIntegrator.hpp
 * @brief Interface for a generic API for mixed integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_MIXEDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_MIXEDINTEGRATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
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
         ~MixedIntegrator() = default;

         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;

         /**
          * @brief Apply FFT
          *
          * @param phys FFT input
          * @param mods FFT output
          */
         void applyFft(Eigen::Ref<MatrixZ> mods, const Eigen::Ref<const Matrix>& phys) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void output(Eigen::Ref<MatrixZ> out) const;

         /**
          * @brief Apply derivative to ouput of given order
          */
         void outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale) const;

         /**
          * @brief Apply derivative to ouput of given order, modified values can be specified for mixed operators
          */
         void outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const;

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
