/**
 * @file AnnulusChebyshevTransform.hpp
 * @brief Implementation of the FFTW transform for a Chebyshev expansion for an annulus radius
 */

#ifndef QUICC_TRANSFORM_ANNULUSCHEBYSHEVTRANSFORM_HPP
#define QUICC_TRANSFORM_ANNULUSCHEBYSHEVTRANSFORM_HPP

// Debug includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Transform/ITransform.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Transform.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the FFTW transform for a Chebyshev expansion for an annulus radius
    */
   class AnnulusChebyshevTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Fft::Chebyshev::LinearMap::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Fft::Chebyshev::LinearMap::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Very basic constructor
          */
         AnnulusChebyshevTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~AnnulusChebyshevTransform();

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
          * @brief Compute forward transform
          *
          * @param rChebVal   Output Chebyshev coefficients
          * @param physVal    Input physical values
          * @param id         Integrator to use
          */
         void forward(MatrixZ& rChebVal, const MatrixZ& physVal, const std::size_t id) override;

         /**
          * @brief Compute backward transform
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          *
          * @param rPhysVal   Output physical values
          * @param chebVal    Input Chebyshev coefficients
          * @param id         Projector to use
          */
         void backward(MatrixZ& rPhysVal, const MatrixZ& chebVal, const std::size_t id) override;

         /**
          * @brief Compute energy reduction operation
          *
          * Compute energy reduction operation
          *
          * @param spectrum   Output energy spectrum
          * @param in         Input spectral coefficients
          * @param id         Energy reductor to use
          */
         void reduce(Matrix& spectrum, const MatrixZ& in, const std::size_t id) override;

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
         Fft::Chebyshev::LinearMap::Transform mImpl;

         //
         // Disabled transforms
         //

         /**
          * @brief Compute forward transform (disabled)
          */
         virtual void forward(Matrix& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute forward transform (disabled)
          */
         virtual void forward(MatrixZ& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute forward transform (disabled)
          */
         virtual void forward(Matrix& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute backward transform (disabled)
          */
         virtual void backward(Matrix& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute backward transform (disabled)
          */
         virtual void backward(MatrixZ& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute backward transform (disabled)
          */
         virtual void backward(Matrix& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute reduction transform (disabled)
          */
         virtual void reduce(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute reduction transform (disabled)
          */
         virtual void reduce(MatrixZ& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute reduction transform (disabled)
          */
         virtual void reduce(Matrix& rOut, const Matrix& in, const std::size_t id) override;
   };

}
}

#endif // QUICC_TRANSFORM_ANNULUSCHEBYSHEVTRANSFORM_HPP
