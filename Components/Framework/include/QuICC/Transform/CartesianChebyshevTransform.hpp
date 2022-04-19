/**
 * @file CartesianChebyshevTransform.hpp
 * @brief Implementation of the FFTW transform for a Chebyshev expansion
 */

#ifndef QUICC_TRANSFORM_CARTESIANCHEBYSHEVTRANSFORM_HPP
#define QUICC_TRANSFORM_CARTESIANCHEBYSHEVTRANSFORM_HPP

// Debug includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//
#include <set>
#include <map>

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
    * @brief Implementation of the FFTW transform for a Chebyshev expansion
    */
   class CartesianChebyshevTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Fft::Chebyshev::LinearMap::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Fft::Chebyshev::LinearMap::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Very basic constructor
          */
         CartesianChebyshevTransform();

         /**
          * @brief Destroy the FFTW plans
          */
         ~CartesianChebyshevTransform();

         /**
          * @brief Initialise the FFT computations (plans, etc)
          *
          * Compute the optimal plans for the required FFT transforms
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

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
          * @brief Compute forward FFT (R2R)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          */
         void forward(Matrix& rOut, const Matrix& in, const std::size_t integrator) override;

         /**
          * @brief Compute forward transform (C2C)
          *
          * Compute the FFT from real physical space to Chebyshev spectral space
          */
         void forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t integrator) override;

         /**
          * @brief Compute backward transform (R2R)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          */
         void backward(Matrix& rOut, const Matrix& in, const std::size_t projector) override;

         /**
          * @brief Compute backward transform (C2C)
          *
          * Compute the FFT from Chebyshev spectral space to real physical space
          */
         void backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t projector) override;

         /**
          * @brief Compute energy reduction operation
          *
          * Compute energy reduction operation
          *
          * @param spectrum   Output energy spectrum
          * @param in         Input spectral coefficients
          * @param reductor   Energy reductor to use
          */
         void reduce(Matrix& spectrum, const MatrixZ& in, const std::size_t reductor) override;

         /**
          * @brief Compute energy reduction operation
          *
          * Compute energy reduction operation
          *
          * @param spectrum   Output energy spectrum
          * @param in         Input spectral coefficients
          * @param reductor   Energy reductor to use
          */
         void reduce(Matrix& spectrum, const Matrix& in, const std::size_t reductor) override;

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
         virtual void forward(MatrixZ& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute forward transform (disabled)
          */
         virtual void forward(Matrix& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute backward transform (disabled)
          */
         virtual void backward(Matrix& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute backward transform (disabled)
          */
         virtual void backward(MatrixZ& rOut, const Matrix& in, const std::size_t id) override;

         /**
          * @brief Compute reduction transform (disabled)
          */
         virtual void reduce(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute reduction transform (disabled)
          */
         virtual void reduce(MatrixZ& rOut, const Matrix& in, const std::size_t id) override;
   };

}
}

#endif // QUICC_TRANSFORM_CARTESIANCHEBYSHEVTRANSFORM_HPP
