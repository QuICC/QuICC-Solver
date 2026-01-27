/**
 * @file IALegendreProjector.hpp
 * @brief Interface for a generic ALegendre FFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_ALEGENDRE_PROJECTOR_IALEGENDREPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_ALEGENDRE_PROJECTOR_IALEGENDREPROJECTOR_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/ALegendre/IALegendreOperator.hpp"
#include "QuICC/Transform/Fft/Backend/ALegendreProjector.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Projector {

   /**
    * @brief Interface for a generic ALegendre FFT based projector
    */
   class IALegendreProjector: public IALegendreOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IALegendreProjector();

         /**
          * @brief Destructor
          */
         virtual ~IALegendreProjector();

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const Matrix& in) const override;

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
         virtual void computeALegendreExpansion(const bool isEven) const = 0;

         /**
          * @brief FFT backend
          */
         Backend::ALegendreProjector mBackend;

      private:
         /**
          * @brief Compute transform  block
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transformBlock(MatrixZ& rOut, const MatrixZ& in, const bool isEven, const bool useReal) const;

         /**
          * @brief Compute transform block
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transformBlock(Matrix& rOut, const Matrix& in, const bool isEven) const;

         /**
          * @brief Apply pre FFT operator
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyPreOperator(const Matrix& in, const bool isEven) const = 0;

         /**
          * @brief Apply post FFT operator
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(Matrix& rOut, const bool isEven) const = 0;

         /**
          * @brief Apply pre FFT operator for component wise operations
          *
          * @param in   Input values
          */
         virtual void applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const = 0;

         /**
          * @brief Apply post FFT operator for component wise operations
          *
          * @param rOut Output values
          */
         virtual void applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const = 0;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(MatrixZ& rOut, const Matrix& in) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP
