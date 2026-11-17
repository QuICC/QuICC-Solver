/**
 * @file ComplexFourierTransform.hpp
 * @brief Implementation of the FFTW transform
 */

#ifndef QUICC_TRANSFORM_COMPLEXFOURIERTRANSFORM_HPP
#define QUICC_TRANSFORM_COMPLEXFOURIERTRANSFORM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Transform/ITransform.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Transform.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Interface class to the FFTW routines
    */
   class ComplexFourierTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Fft::Fourier::Complex::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Fft::Fourier::Complex::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Very basic constructor
          */
         ComplexFourierTransform() = default;

         /**
          * @brief Destroy the FFTW plans
          */
         ~ComplexFourierTransform() = default;

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
          * @brief Compute forward FFT (C2C)
          *
          * Compute the FFT from complex physical space to complex spectral space
          */
         void forward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id) override;

         /**
          * @brief Compute backward FFT (C2C)
          *
          * Compute the FFT from complex spectral space to complex physical space
          */
         void backward(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const std::size_t id) override;

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
         Fft::Fourier::Complex::Transform mImpl;
   };

}
}

#endif // QUICC_TRANSFORM_COMPLEXFOURIERTRANSFORM_HPP
