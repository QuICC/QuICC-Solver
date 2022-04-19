/** 
 * @file IWorlandEnergy.hpp
 * @brief Interface for a generic Worland FFT based energy operator 
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_IWORLANDENERGY_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_IWORLANDENERGY_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/IWorlandOperator.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/WorlandProjector.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/WorlandIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a generic Worland FFT based energy operator
    */ 
   class IWorlandEnergy: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandEnergy(const int shift);

         /**
          * @brief Destructor
          */
         virtual ~IWorlandEnergy();

         /**
          * @brief Compute transform R2R componentwise
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void transform(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform R2R
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
          * @brief Compute energy quadrature
          */
         void computeEnergyQuadrature(internal::Array& igrid, internal::Array& iweights, Array& eweights, const int iL) const;

         /**
          * @brief Initialise FFT backend
          */
         virtual void initBackend() const override;

         /**
          * @brief Polynomial shift
          */
         const int mcShift;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable internal::Array  mWeights;

         /**
          * @brief Storage for the energy weights
          */
         mutable internal::Array  mEWeights;

         /**
          * @brief Storage for the projector
          */
         mutable std::vector<Matrix>  mPOps;

         /**
          * @brief Storage for the energy integrator
          */
         mutable std::vector<Matrix>  mEOps;

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const Array& eweights, const int i) const = 0;

         /**
          * @brief Compute energy (integral of squared values)
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         void applyOperators(Matrix& rOut, const MatrixZ& in) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const = 0;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_IWORLANDENERGY_HPP
