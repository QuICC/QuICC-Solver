/**
 * @file IWorlandProjector.hpp
 * @brief Interface for a Worland based projector
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP

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
#include "QuICC/Transform/Poly/Worland/IWorlandOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Projector {

   /**
    * @brief Interface for a Worland based projector
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
         virtual ~IWorlandProjector();

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
          * @brief Default implementation to apply ith operator
          */
         void defaultApplyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;

         /**
          * @brief Storage for the operators
          */
         mutable std::vector<Matrix>  mOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable internal::Array  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable internal::Array  mWeights;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const internal::Array& igrid, const internal::Array& iweights) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const internal::Array& igrid, const internal::Array& iweights, const int i) const = 0;

         /**
          * @brief Compute polynomial projection
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;

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
         virtual void applyOperator(Eigen::Ref<MatrixZ> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const = 0;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_WORLAND_PROJECTOR_IWORLANDPROJECTOR_HPP
