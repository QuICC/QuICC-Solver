/**
 * @file ILinearMapProjector.hpp
 * @brief Interface for a generic Chebyshev FFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ILINEARMAPPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ILINEARMAPPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Backend/ChebyshevProjector.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

/**
 * @brief Interface for a generic Chebyshev FFT based projector
 */
class ILinearMapProjector : public IChebyshevOperator
{
public:
   /**
    * @brief Constructor
    */
   ILinearMapProjector() = default;

   /**
    * @brief Destructor
    */
   virtual ~ILinearMapProjector() = default;

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
    * @brief FFT backend
    */
   Backend::ChebyshevProjector mBackend;

private:
   /**
    * @brief Apply pre FFT operator
    *
    * @param tmp  Temporary padded modal values
    * @param in   Input values
    */
   virtual void applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const = 0;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   virtual void applyPostOperator(Eigen::Ref<Matrix> rOut) const = 0;

   /**
    * @brief Apply pre FFT operator for component wise operations
    *
    * @param tmp Temporary padded modal values, either real or im part only
    * @param in Input values
    * @param useReal 1 -> extract real part, 0 -> extract im part
    */
   virtual void applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
      const bool useReal) const = 0;

   /**
    * @brief Apply post FFT operator for component wise operations
    *
    * @param rOut Output values
    */
   virtual void applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
      const bool useReal) const = 0;

   /**
    * @brief Compute transform R2C (disabled)
    *
    * @param rOut Output values
    * @param in   Input values
    */
   virtual void transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in) const override;

   /**
    * @brief Compute transform C2R (disabled)
    *
    * @param rOut Output values
    * @param in   Input values
    */
   virtual void transform(Eigen::Ref<Matrix> rOut, const Eigen::Ref<const MatrixZ>& in) const override;
};

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ILINEARMAPPROJECTOR_HPP
