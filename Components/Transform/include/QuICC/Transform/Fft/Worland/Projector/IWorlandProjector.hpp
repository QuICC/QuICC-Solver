/**
 * @file IWorlandProjector.hpp
 * @brief Interface for a generic Worland FFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/IWorlandOperator.hpp"
#include "QuICC/Transform/Fft/Backend/WorlandProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   /**
    * @brief Interface for a generic Worland FFT based projector
    */
   class IWorlandProjector: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandProjector();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandProjector() = default;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const MatrixZ>& in) const override;

         /**
          * @brief Compute transform
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
          * @brief FFT backend
          */
         Backend::WorlandProjector mBackend;

      private:
         /**
          * @brief Compute transform  block
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

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP
