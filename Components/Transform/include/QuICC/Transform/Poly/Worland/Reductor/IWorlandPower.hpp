/**
 * @file IWorlandPower.hpp
 * @brief Interface for a Worland based power operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDPOWER_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDPOWER_HPP

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
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandReductor.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   /**
    * @brief Interface for a Worland based power operator
    */
   class IWorlandPower: public IWorlandReductor
   {
      public:
         /**
          * @brief Constructor
          */
         IWorlandPower(const int shift);

         /**
          * @brief Destructor
          */
         virtual ~IWorlandPower();

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
         void computePowerQuadrature(internal::Array& igrid, internal::Array& iweights, const int gSize) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const = 0;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const internal::Array& igrid, const internal::Array& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const int i) const = 0;

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

#endif // QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_IWORLANDPOWER_HPP
