/**
 * @file IOperator.hpp
 * @brief Interface for a Finite Differences sphere operator
 */

#ifndef QUICC_TRANSFORM_FINITEDIFF_SPHERE_IOPERATOR_HPP
#define QUICC_TRANSFORM_FINITEDIFF_SPHERE_IOPERATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Transform/FiniteDiff/Setup.hpp"
#include "QuICC/Transform/ITransformOperator.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

   /**
    * @brief Interface for a Finite Differences sphere operator
    */
   class IOperator: public ITransformOperator
   {
      public:
         /// Typedef for the configuration class
         typedef Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         IOperator();

         /**
          * @brief Destructor
          */
         virtual ~IOperator() = default;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup, const Internal::Array& igrid) const;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const override;

         /**
          * @brief Initialise the transform
          *
          * @param spSetup   Shared setup object for the transform
          */
         void init(SharedTransformSetup spSetup) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(MatrixZ& rOut, const MatrixZ& in) const override;

         /**
          * @brief Compute transform
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void transform(Matrix& rOut, const MatrixZ& in) const;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const override;

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const = 0;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const = 0;

      protected:
         /**
          * @brief Check grid size can accommodate exact computation
          */
         void checkGridSize(const int n, const int l, const int nG) const;

         /**
          * @brief Polynomial setup object providing the sizes
          */
         mutable SharedSetup    mspSetup;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const Internal::Array& igrid) const = 0;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(MatrixZ& rOut, const MatrixZ& in) const;

         /**
          * @brief Apply operators
          *
          * @param rOut Output values
          * @param in   Input values
          */
         virtual void applyOperators(Matrix& rOut, const MatrixZ& in) const;
   };

} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FINITEDIFF_SPHERE_IOPERATOR_HPP
