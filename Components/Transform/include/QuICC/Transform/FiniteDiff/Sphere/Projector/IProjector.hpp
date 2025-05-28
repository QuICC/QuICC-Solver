/**
 * @file IProjector.hpp
 * @brief Interface for a Finite Differences Sphere projector
 */

#ifndef QUICC_TRANSFORM_FINITEDIFF_SPHERE_PROJECTOR_IPROJECTOR_HPP
#define QUICC_TRANSFORM_FINITEDIFF_SPHERE_PROJECTOR_IPROJECTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Projector {

   /**
    * @brief Interface for a Finite Differences Sphere projector
    */
   class IProjector: public IOperator
   {
      public:

         using OpMatrixR = Eigen::Ref<MatrixZ>;;
         using OpMatrixCR = Eigen::Ref<const MatrixZ>;

         /**
          * @brief Constructor
          */
         IProjector();

         /**
          * @brief Destructor
          */
         virtual ~IProjector() = default;

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
         void defaultApplyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

         /**
          * @brief Storage for the operators
          */
         mutable std::vector<Matrix>  mOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable Internal::Array  mGrid;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const Internal::Array& igrid) const override;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, const Internal::Array& igrid, const int i) const = 0;

         /**
          * @brief Compute backward transform
          *
          * @param rOut Output physical values
          * @param in   Input spectral coefficients
          */
         void applyOperators(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute backward transform
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

#endif // QUICC_TRANSFORM_FINITEDIFF_SPHERE_PROJECTOR_IPROJECTOR_HPP
