/**
 * @file SphereFftWorlandTransform.hpp
 * @brief Implementation of the FFT based Worland transform in a sphere
 */

#ifndef QUICC_TRANSFORM_SPHEREFFTWORLANDTRANSFORM_HPP
#define QUICC_TRANSFORM_SPHEREFFTWORLANDTRANSFORM_HPP

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
#include "QuICC/Transform/Fft/Worland/Transform.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of the Worland transform in a sphere
    */
   class SphereFftWorlandTransform: public ITransform
   {
      public:
         /// Typedef for setup type
         typedef Fft::Worland::Transform::SetupType SetupType;

         /// Typedef for setup type
         typedef Fft::Worland::Transform::SharedSetupType SharedSetupType;

         /**
          * @brief Constructor
          */
         SphereFftWorlandTransform();

         /**
          * @brief Destructor
          */
         virtual ~SphereFftWorlandTransform();

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
          * @brief Initialise the polynomial transform (matrices, weights, grid, etc)
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedSetupType spSetup);

         /**
          * @brief Compute quadrature integration
          *
          * @param rOut       Output spectral coefficients
          * @param in         Input physical values
          * @param id         Integrator to use
          */
         virtual void forward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute polynomial projection
          *
          * @param rOut       Output physical values
          * @param in         Input spectral coefficients
          * @param id         Projector to use
          */
         virtual void backward(MatrixZ& rOut, const MatrixZ& in, const std::size_t id) override;

         /**
          * @brief Compute reduction operation
          *
          * @param spectrum   Output energy spectrum
          * @param in         Input spectral coefficients
          * @param id         Energy reductor to use
          */
         virtual void reduce(Matrix& spectrum, const MatrixZ& in, const std::size_t id) override;

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
         Fft::Worland::Transform mImpl;

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

#endif // QUICC_TRANSFORM_SPHEREFFTWORLANDTRANSFORM_HPP
