/**
 * @file IBesselPower.hpp
 * @brief Interface for a Bessel based power operator
 */

#ifndef QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_IBESSELPOWER_HPP
#define QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_IBESSELPOWER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/IBesselReductor.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Bessel {

namespace Reductor {

   /**
    * @brief Interface for a Bessel based power operator
    */
   class IBesselPower: public IBesselReductor
   {
      public:
         /**
          * @brief Constructor
          */
         IBesselPower(const int shift);

         /**
          * @brief Destructor
          */
         virtual ~IBesselPower() = default;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

      protected:
         /**
          * @brief Apply ith operator
          */
         virtual void defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief Storage for the power integrator
          */
         mutable std::vector<Matrix>  mEOps;

         /**
          * @brief Polynomial shift
          */
         const int mcShift;

         /**
          * @brief Compute power quadrature
          */
         void computePowerQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int gSize) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const = 0;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const = 0;

         /**
          * @brief Compute power (integral of squared values)
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         virtual void applyOperators(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute power (integral of squared values)
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_BESSEL_REDUCTOR_IBESSELPOWER_HPP
