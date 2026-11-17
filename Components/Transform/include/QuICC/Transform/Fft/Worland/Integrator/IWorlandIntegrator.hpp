/**
 * @file IWorlandIntegrator.hpp
 * @brief Interface for a generic Worland FFT based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_IWORLANDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_IWORLANDINTEGRATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/IWorlandOperator.hpp"
#include "QuICC/Transform/Fft/Backend/WorlandIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   /**
    * @brief Interface for a generic Worland FFT based integrator
    */
   class IWorlandIntegrator: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandIntegrator() = default;

         /**
          * @brief Compute transform R2R componentwise
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Compute transform R2R
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in) const override;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

      protected:
         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief Initialise FFT backend
          */
         virtual void computeWorlandExpansion(const bool isEven) const = 0;

         /**
          * @brief Compute transform block
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transformBlock(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in, const bool isEven, const bool useReal) const;

         /**
          * @brief Compute transform block
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transformBlock(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const Matrix>& in, const bool isEven) const;

         /**
          * @brief FFT backend
          */
         Backend::WorlandIntegrator mBackend;

      private:
         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(const Eigen::Ref<const Matrix>& in, const bool isEven) const = 0;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Eigen::Ref<Matrix> rOut, const bool isEven) const = 0;

         /**
          * @brief Apply pre FFT operator for component wise operations
          *
          * @param in   Input values
          */
         virtual void applyPreOperator(const Eigen::Ref<const MatrixZ>& in, const bool isEven, const bool useReal) const = 0;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Eigen::Ref<MatrixZ> rOut, const bool isEven, const bool useReal) const = 0;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_INTEGRATOR_IWORLANDINTEGRATOR_HPP
