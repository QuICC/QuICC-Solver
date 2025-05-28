/**
 * @file IIntegrator.hpp
 * @brief Interface for a Finite Difference sphere integrator
 */

#ifndef QUICC_TRANSFORM_FINITEDIFF_SPHERE_INTEGRATOR_IINTEGRATOR_HPP
#define QUICC_TRANSFORM_FINITEDIFF_SPHERE_INTEGRATOR_IINTEGRATOR_HPP

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

namespace Integrator {

   /**
    * @brief Interface for a Finite Difference sphere integrator
    */
   class IIntegrator: public IOperator
   {
      public:

         using OpMatrixR = Eigen::Ref<MatrixZ>;
         using OpMatrixCR = Eigen::Ref<const MatrixZ>;

         /**
          * @brief Constructor
          */
         IIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~IIntegrator() = default;

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
          * @brief Compute forward transform
          *
          * @param rOut Output spectral coefficients
          * @param in   Input physical values
          */
         void applyOperators(MatrixZ& rPhysVal, const MatrixZ& specVal) const override;

         /**
          * @brief Compute forward transform
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

} // namespace Integrator
} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FINITEDIFF_SPHERE_IINTEGRATOR_HPP
