/**
 * @file ILinearMapEnergy.hpp
 * @brief Interface for a generic Chebyshev FFT based energy reductor
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ILINEARMAPENERGY_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ILINEARMAPENERGY_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Backend/ChebyshevEnergy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"
#include "Types/Typedefs.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

/**
 * @brief Interface for a generic Chebyshev FFT based energy reductor
 */
class ILinearMapEnergy : public IChebyshevOperator
{
public:
   /**
    * @brief Constructor
    */
   ILinearMapEnergy() = default;

   /**
    * @brief Destructor
    */
   virtual ~ILinearMapEnergy() = default;

   /**
    * @brief Compute reduction of complex data
    *
    * @param rOut Output values
    * @param in   Input values
    */
   virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

   /**
    * @brief Compute reduction of real data, anelastic case
    *
    * @param rOut Output values
    * @param in   Input values
    * @param pF   Radial profile (e.g density)
    */
   virtual void transform(Matrix& rOut, const MatrixZ& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const override;


   /**
    * @brief Compute reduction of real data
    *
    * @param rOut Output values
    * @param in   Input values
    */
   virtual void transform(Matrix& rOut, const Matrix& in) const override;

    /**
    * @brief Compute reduction of complex data, anelastic case
    *
    * @param rOut Output values
    * @param in   Input values
    * @param pF   Radial profile (e.g density)
    */
   virtual void transform(Matrix& rOut, const Matrix& in, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const override;

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

   /**
    * @brief Unused interface
    */
   void transform(MatrixZ& rOut, const MatrixZ& in) const final;

protected:
   /**
    * @brief Initialise FFT backend
    */
   virtual void initBackend() const override;

   /**
    * @brief Initialise FFT backend, anelastic overload
    * @param pF   Shared pointer to radial profile
    */
   virtual void initBackendAnelastic(std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const override;


   /**
    * @brief FFT backend
    */
   Backend::ChebyshevEnergy mBackend;

   mutable int mExtraSize = 0;

private:
   /**
    * @brief Apply pre FFT operator
    *
    * @param rOut Output values
    * @param in   Input values
    */
   virtual void applyPreOperator(Matrix& tmp, const Matrix& in) const = 0;

   /**
    * @brief Apply post FFT operator
    *
    * @param rOut Output values
    */
   virtual void applyPostOperator(Matrix& rOut, const Matrix& tmp) const = 0;

   /**
    * @brief Apply pre FFT operator for component wise operations
    *
    * @param in   Input values
    */
   virtual void applyPreOperator(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const = 0;
};

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_LINEARMAP_ILINEARMAPENERGY_HPP
