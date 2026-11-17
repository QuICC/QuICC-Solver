/**
 * @file DivY1D1Y1.hpp
 * @brief Implementation of the Chebyshev based 1/Y D Y projector, with linear
 * map y = ax + b
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D1Y1_BASE_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D1Y1_BASE_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/ILinearMapProjector.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Tags.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

template <typename Impl> class DivY1D1Y1;

/**
 * @brief Implementation of the Chebyshev based 1/Y D Y projector, with linear
 * map y = ax +  b
 */
template <> class DivY1D1Y1<base_t> : public ILinearMapProjector
{
public:
   /**
    * @brief Constructor
    */
   DivY1D1Y1() = default;

   /**
    * @brief Destructor
    */
   ~DivY1D1Y1() = default;

protected:
   /**
    * @brief Initialize storage
    */
   void initBackend() const final;

private:
   /**
    * @brief Initialize solver operators
    */
   void initOperator() const final;

   /**
    * @brief Apply pre FFT operator
    *
    * @param tmp Temporary padded modal values
    * @param in   Input values
    */
   void applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const final;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   void applyPostOperator(Eigen::Ref<Matrix> rOut) const final;

   /**
    * @brief Apply pre FFT operator for component wise openerations
    *
    * @param tmp Temporary padded modal values, either real or im part only
    * @param in Input values
    * @param useReal Real vs Imag flag
    */
   void applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
      const bool useReal) const final;

   /**
    * @brief Apply post FFT operator for component wise operations
    *
    * @param rOut Output values
    * @param useReal Real vs Imag flag
    */
   void applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
      const bool useReal) const final;
};

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_PROJECTOR_DIVY1D1Y1_BASE_HPP
