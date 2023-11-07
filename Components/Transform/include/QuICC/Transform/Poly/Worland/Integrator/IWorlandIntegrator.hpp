/**
 * @file IWorlandIntegrator.hpp
 * @brief Interface for a Worland based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_IWORLANDINTEGRATOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_INTEGRATOR_IWORLANDINTEGRATOR_HPP


// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Integrator {

   /**
    * @brief Interface for a Worland based integrator
    */
   class IWorlandIntegrator: public IWorlandOperator
   {
      public:

         using OpMatrixR = Eigen::Ref<MatrixZ>;
         using OpMatrixCR = Eigen::Ref<const MatrixZ>;

         /**
          * @brief Constructor
          */
         IWorlandIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandIntegrator();

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

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
          * @brief Truncate Quasi-inverse ?
          */
         const bool mcTruncQI;

         /**
          * @brief Default implementation to apply ith operator
          */
         void defaultApplyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

         /**
          * @brief Storage for the operators
          */
         mutable std::vector<Matrix>  mOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable Internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable Internal::Array  mWeights;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const Internal::Array& iweights, const int i) const = 0;
         /**
          * @brief Compute polynomial integration
          *
          * @param rOut Output spectral coefficients
          * @param in   Input physical values
          */
         void applyOperators(MatrixZ& rPhysVal, const MatrixZ& specVal) const override;

         /**
          * @brief Compute polynomial projection
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         void applyOperators(Matrix& rOut, const MatrixZ& in) const override;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const = 0;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_IWORLANDINTEGRATOR_HPP
