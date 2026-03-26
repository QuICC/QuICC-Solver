/**
 * @file IChebyshevOperator.hpp
 * @brief Interface for a generic Chebyshev FFT based operator
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVOPERATOR_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Chebyshev/Setup.hpp"
#include "QuICC/Transform/ITransformOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

   /**
    * @brief Interface for a generic Chebyshev FFT based operator
    */
   class IChebyshevOperator: public ITransformOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IChebyshevOperator() = default;

         /**
          * @brief Destructor
          */
         virtual ~IChebyshevOperator() = default;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup) const override;

         /**
          * @brief Initialise the transform, anelastic overload
          *
          * @param spSetup   Shared setup object for the transform
          * @param pF        Shared pointer to radial profile
          */
         void init(SharedTransformSetup spSetup, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const override;

         using ITransformOperator::transform;

         /**
          * @brief Compute Chebyshev transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const = 0;

         /**
          * @brief Compute transform C2R, anelastic version
          *
          * @param rOut Output values
          * @param in   Input values
          * @param pF   radial profile
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
         {
            // Default implementation just ignores pF and calls the standard transform
            this->transform(rOut, in);
         }

         /**
          * @brief Compute Chebyshev transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const = 0;

         /**
          * @brief Compute transform R2C, anelastic version
          *
          * @param rOut Output values
          * @param in   Input values
          * @param pF   radial profile
          */
         virtual void transform(Matrix& rOut, const Matrix& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
         {
            // Default implementation just ignores pF and calls the standard transform
            this->transform(rOut, in);
         }

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const = 0;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const = 0;

         /**
          * @brief Cleanup
          */
         void cleanup();

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

      protected:
         /**
          * @brief Initialize base
          */
         void initBase() const;

         /**
          * @brief Initialize base, anelastic overload
          */
         void initBase(std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

      private:
         /**
          * @brief Operator specific initialization with empty default implementation
          */
         virtual void initOperator() const;

         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const = 0;

         /**
          * @brief Initialise FFT backend, anelastic overload
          * @param pF   Shared pointer to radial profile
          */
         virtual void initBackendAnelastic(std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
         {
            // Default implementation: just call standard initBackend
            this->initBackend();
         }
   };

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_ICHEBYSHEVOPERATOR_HPP
