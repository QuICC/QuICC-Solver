/**
 * @file MixedFourierTransform.hpp
 * @brief Implementation of the FFTW mixed transform
 */

#ifndef QUICC_TRANSFORM_MIXEDFOURIERTRANSFORM_HPP
#define QUICC_TRANSFORM_MIXEDFOURIERTRANSFORM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Transform/ITransform.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Transform.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Interface class to the FFTW routines
    */
   class MixedFourierTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Fft::Fourier::Mixed::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Fft::Fourier::Mixed::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Ctor
          */
         MixedFourierTransform() = default;

         /**
          * @brief Dtor
          */
         ~MixedFourierTransform() = default;

         /**
          * @brief set list of required options
          */
         virtual void requiredOptions(std::set<std::size_t>& list, const Dimensions::Transform::Id dimId) const override;

         /**
          * @brief Set the required options
          */
         virtual void setOptions(const std::map<std::size_t, NonDimensional::SharedINumber>& options, const Dimensions::Transform::Id dimId) override;

         /**
          * @brief Get the physical grid
          */
         virtual Array meshGrid() const override;

         /**
          * @brief Initialise the FFT computations (plans, etc)
          *
          * Compute the optimal plans for the required FFT transforms
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Compute forward FFT (R2C)
          *
          * Compute the FFT from real physical space to complex spectral space
          */
         void forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in, const std::size_t id) override;

         /**
          * @brief Compute backward FFT (C2R)
          *
          * Compute the FFT from complex spectral space to real physical space
          */
         void backward(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id) override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

         /**
          * @brief Profile the memory requirements
          */
         virtual void profileStorage() const override;

      protected:

      private:
         /**
          * @brief Initialise the operators
          */
         void initOperators();

         /**
          * @brief Transform implementation
          */
         Fft::Fourier::Mixed::Transform mImpl;
   };

}
}

#endif // QUICC_TRANSFORM_MIXEDFOURIERTRANSFORM_HPP
